/*
 * $Id: OperatorImageTileStack.java,v 1.1 2010-01-14 17:22:47 lveci Exp $
 *
 * Copyright (C) 2010 by Brockmann Consult (info@brockmann-consult.de)
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation. This program is distributed in the hope it will
 * be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */
package org.esa.beam.framework.gpf.internal;

import com.bc.ceres.core.Assert;

import org.esa.beam.framework.datamodel.Band;
import org.esa.beam.framework.gpf.Tile;
import org.esa.beam.jai.ImageManager;
import org.esa.beam.util.ImageUtils;

import java.awt.Dimension;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.image.Raster;
import java.awt.image.SampleModel;
import java.awt.image.WritableRaster;
import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;

import javax.media.jai.PlanarImage;


class OperatorImageTileStack extends OperatorImage {

    private final Object[][] locks;
    
    OperatorImageTileStack(Band targetBand, OperatorContext operatorContext, Object[][] locks) {
        super(targetBand, operatorContext);
        this.locks = locks;
    }
    
    @Override
    public Raster computeTile(int tileX, int tileY) {
        Object tileLock = locks[tileX][tileY];
        
        // lock to prevent multiple simultaneous computations.
        synchronized (tileLock) {
            Raster tileFromCache = getTileFromCache(tileX, tileY);
            if (tileFromCache != null) {
                return tileFromCache;
            } else {
                /* Create a new WritableRaster to represent this tile. */
                Point location = new Point(tileXToX(tileX), tileYToY(tileY));
                WritableRaster dest = createWritableRaster(sampleModel, location);
                
                /* Clip output rectangle to image bounds. */
                Rectangle rect = new Rectangle(location.x, location.y,
                                               sampleModel.getWidth(),
                                               sampleModel.getHeight());
                Rectangle destRect = rect.intersection(getBounds());
                computeRect((PlanarImage[])null, dest, destRect);
                return dest;
            }
        }
    }

    @Override
    protected void computeRect(PlanarImage[] ignored, WritableRaster tile, Rectangle destRect) {

        long nanos1 = System.nanoTime();

        OperatorContext operatorContext = getOperatorContext();
        Band[] targetBands = operatorContext.getTargetProduct().getBands();
        Map<Band, Tile> targetTiles = new HashMap<Band, Tile>(targetBands.length * 2);
        if (operatorContext.isPassThrough()) {
            for (Band band : targetBands) {
                if (isBandComputedByThisOperator(band)) {
                    targetTiles.put(band, operatorContext.getSourceTile(band, destRect, getProgressMonitor()));
                }
            }
            operatorContext.getOperator().computeTileStack(targetTiles, destRect, getProgressMonitor());
        } else {
            Map<Band, WritableRaster> writableRasters = new HashMap<Band, WritableRaster>(targetBands.length);
            for (Band band : targetBands) {
                if (isBandComputedByThisOperator(band)) {
                    WritableRaster tileRaster = getWritableRaster(band, tile);
                    writableRasters.put(band, tileRaster);
                    targetTiles.put(band, createTargetTile(band, tileRaster, destRect));
                }
            }
            operatorContext.getOperator().computeTileStack(targetTiles, destRect, getProgressMonitor());
            
            for (Entry<Band, WritableRaster> entry : writableRasters.entrySet()) {
                Band band = entry.getKey();
                WritableRaster writableRaster = entry.getValue();
                // casting to access "addTileToCache" method
                OperatorImageTileStack operatorImage = (OperatorImageTileStack) operatorContext.getTargetImage(band);
                final int tileX = XToTileX(destRect.x);
                final int tileY = YToTileY(destRect.y);
                //put raster into cache after computing them.
                operatorImage.addTileToCache(tileX, tileY, writableRaster);
            }
        }
        long nanos2 = System.nanoTime();
        updatePerformanceMetrics(nanos1, nanos2, destRect);
    }
    
    private WritableRaster getWritableRaster(Band band, WritableRaster targetTileRaster) {
        WritableRaster tileRaster;
        if (band == getTargetBand()) {
            tileRaster = targetTileRaster;
        } else {
            OperatorContext operatorContext = getOperatorContext();
            // casting to access "getWritableRaster" method
            OperatorImageTileStack operatorImage = (OperatorImageTileStack) operatorContext.getTargetImage(band);
            Assert.state(operatorImage != this);
            tileRaster = operatorImage.getWritableRaster(targetTileRaster.getBounds());
        }
        return tileRaster;
    }

    private WritableRaster getWritableRaster(Rectangle tileRectangle) {
        Assert.argument(tileRectangle.x % getTileWidth() == 0, "rectangle");
        Assert.argument(tileRectangle.y % getTileHeight() == 0, "rectangle");
        Assert.argument(tileRectangle.width == getTileWidth(), "rectangle");
        Assert.argument(tileRectangle.height == getTileHeight(), "rectangle");
        final int tileX = XToTileX(tileRectangle.x);
        final int tileY = YToTileY(tileRectangle.y);
        Raster tileFromCache = getTileFromCache(tileX, tileY);
        WritableRaster writableRaster;
        if (tileFromCache != null) {
            // we already have a WritableRaster in the cache
            writableRaster = (WritableRaster) tileFromCache;
        } else {
            writableRaster = createWritableRaster(tileRectangle);
        }
        return writableRaster;
    }
    
    private WritableRaster createWritableRaster(Rectangle rectangle) {
        final int dataBufferType = ImageManager.getDataBufferType(getTargetBand().getDataType());
        SampleModel sampleModel = ImageUtils.createSingleBandedSampleModel(dataBufferType, rectangle.width,
                                                                           rectangle.height);
        final Point location = new Point(rectangle.x, rectangle.y);
        return createWritableRaster(sampleModel, location);
    }
    
    private boolean isBandComputedByThisOperator(Band band) {
        if (band == getTargetBand()) {
            return true;
        }
        if (!band.isSourceImageSet()) {
            return false;
        }
        OperatorImage image = getOperatorContext().getTargetImage(band);
        return image != null && image == band.getSourceImage().getImage(0);
    }

    /**
     * Create a lock objects for each tile. These locks are used by all images in the tile stack.
     * This prevent multiple computation of tiles. 
     */
    static Object[][] createLocks(int width, int height, Dimension tileSize) {
        int tw = tileSize.width;
        int numXTiles =  PlanarImage.XToTileX(0 + width - 1, 0, tw) - PlanarImage.XToTileX(0, 0, tw) + 1;
        int th = tileSize.height;
        int numYTiles =  PlanarImage.YToTileY(0 + height - 1, 0, th) - PlanarImage.YToTileY(0, 0, th) + 1;
        final Object[][] lock = new Object[numXTiles][numYTiles];
        for (int x = 0; x < numXTiles; x++) {
            for (int y = 0; y < numYTiles; y++) {
                lock[x][y] = new Object();
            }
        }
        return lock;
    }
}
