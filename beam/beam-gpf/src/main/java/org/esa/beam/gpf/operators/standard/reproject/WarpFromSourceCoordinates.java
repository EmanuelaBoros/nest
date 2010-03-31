/*
 * $Id: WarpFromSourceCoordinates.java,v 1.1 2010-03-31 13:59:24 lveci Exp $
 *
 * Copyright (C) 2009 by Brockmann Consult (info@brockmann-consult.de)
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
package org.esa.beam.gpf.operators.standard.reproject;

import java.awt.Rectangle;
import java.awt.image.DataBuffer;
import java.awt.image.Raster;

import javax.media.jai.OpImage;
import javax.media.jai.PixelAccessor;
import javax.media.jai.UnpackedImageData;
import javax.media.jai.Warp;

/**
 * Gives the warp coordinates using the data from the given image.
 *
 * @author Marco Zuehlke
 * @since BEAM 4.7
 */
class WarpFromSourceCoordinates extends Warp {
    
    private final OpImage opImage;

    WarpFromSourceCoordinates(OpImage opImage) {
        this.opImage = opImage;
    }

    /**
     * Computes the source pixel positions for a given rectangular
     * destination region, subsampled with an integral period.
     */
    @Override
    public float[] warpSparseRect(int xmin, int ymin, int width, int height, int periodX, int periodY, float[] destRect) {
        if (periodX < 1) throw new IllegalArgumentException(String.valueOf(periodX));
        if (periodY < 1) throw new IllegalArgumentException(String.valueOf(periodY));
        
        final int xmax  = xmin + width;
        final int ymax  = ymin + height;
        final int count = ((width+(periodX-1))/periodX) * ((height+(periodY-1))/periodY);
        if (destRect == null) {
            destRect = new float[2*count];
        }
        
        Rectangle bounds = new Rectangle(xmin, ymin, width, height);
        int xIDNew = opImage.XToTileX(xmin);
        int yIDNew = opImage.YToTileY(ymin);
        Raster tile = opImage.getTile(xIDNew, yIDNew);
        PixelAccessor accessor = new PixelAccessor(opImage);
        UnpackedImageData srcImD = accessor.getPixels(tile, bounds, DataBuffer.TYPE_FLOAT, false);
        float[] data = srcImD.getFloatData(0);
        int sxStart = srcImD.bandOffsets[0];
        int syStart = srcImD.bandOffsets[1];
        int pixelStride = srcImD.pixelStride;
        int lineStride = srcImD.lineStride;

        int index = 0;
        for (int y=ymin; y<ymax; y+=periodY) {
            int sxPos = sxStart;
            int syPos = syStart;
            for (int x=xmin; x<xmax; x+=periodX) {
                destRect[index++] = data[sxPos];
                destRect[index++] = data[syPos];
                sxPos += pixelStride;
                syPos += pixelStride;
            }
            sxStart += lineStride;
            syStart += lineStride;
        }
        return destRect;
    }
}
