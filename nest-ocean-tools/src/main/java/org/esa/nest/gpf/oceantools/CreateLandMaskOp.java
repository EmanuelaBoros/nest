/*
 * Copyright (C) 2010 Array Systems Computing Inc. http://www.array.ca
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation; either version 3 of the License, or (at your option)
 * any later version.
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, see http://www.gnu.org/licenses/
 */
package org.esa.nest.gpf.oceantools;

import com.bc.ceres.core.ProgressMonitor;
import org.esa.beam.framework.datamodel.*;
import org.esa.beam.framework.dataop.dem.ElevationModel;
import org.esa.beam.framework.dataop.dem.ElevationModelDescriptor;
import org.esa.beam.framework.dataop.dem.ElevationModelRegistry;
import org.esa.beam.framework.dataop.resamp.ResamplingFactory;
import org.esa.beam.framework.gpf.Operator;
import org.esa.beam.framework.gpf.OperatorException;
import org.esa.beam.framework.gpf.OperatorSpi;
import org.esa.beam.framework.gpf.Tile;
import org.esa.beam.framework.gpf.annotations.OperatorMetadata;
import org.esa.beam.framework.gpf.annotations.Parameter;
import org.esa.beam.framework.gpf.annotations.SourceProduct;
import org.esa.beam.framework.gpf.annotations.TargetProduct;
import org.esa.beam.util.ProductUtils;
import org.esa.nest.gpf.OperatorUtils;

import java.awt.*;
import java.util.ArrayList;
import java.util.Map;
import java.util.Set;

/**
 * The CreateLandMask operator.
 */
@OperatorMetadata(alias = "Create-LandMask",
        category = "Ocean-Tools",
        description = "Creates a bitmask defining land vs ocean.")
public class CreateLandMaskOp extends Operator {

    @SourceProduct(alias="source")
    private Product sourceProduct;
    @TargetProduct
    private Product targetProduct = null;

    @Parameter(description = "The list of source bands.", alias = "sourceBands", itemAlias = "band",
            sourceProductId = "source", label = "Source Bands")
    private String[] sourceBandNames = null;

    @Parameter(label="Mask the Land", defaultValue = "true")
    private boolean landMask = true;

    @Parameter(label="Use SRTM 3sec", defaultValue = "true")
    private boolean useSRTM = true;

    @Parameter(label="Geometry", defaultValue = "")
    private String geometry = "";

    @Parameter(label="Invert Geometry", defaultValue = "false")
    private boolean invertGeometry = false;

    @Parameter(label="Bypass", defaultValue = "false")
    private boolean byPass = false;

    private ElevationModel dem = null;
    private final static int landThreshold = -10;
    private final static int seaThreshold = -10;

    private final static String tmpVirtBandName = "_tmpVirtBand";

    @Override
    public void initialize() throws OperatorException {
        try {

            targetProduct = new Product(sourceProduct.getName(), sourceProduct.getProductType(),
                    sourceProduct.getSceneRasterWidth(), sourceProduct.getSceneRasterHeight());

            OperatorUtils.copyProductNodes(sourceProduct, targetProduct);

            addSelectedBands();

        } catch (Throwable e) {
            OperatorUtils.catchOperatorException(getId(), e);
        }
    }

    /**
     * Add the user selected bands to target product.
     * @throws OperatorException The exceptions.
     */
    private void addSelectedBands() throws OperatorException {

        final Band[] sourceBands = OperatorUtils.getSourceBands(sourceProduct, sourceBandNames);
        for (Band srcBand : sourceBands) {

            if(geometry != null && !geometry.isEmpty() && !byPass) {
                String expression = geometry+" ? "+srcBand.getName()+" : "+srcBand.getNoDataValue();
                if(invertGeometry) {
                    expression = "!"+expression;
                }
                final VirtualBand virtBand = new VirtualBand(srcBand.getName() + tmpVirtBandName,
                        srcBand.getDataType(),
                        srcBand.getSceneRasterWidth(), srcBand.getSceneRasterHeight(),
                        expression);
                virtBand.setUnit(srcBand.getUnit());
                virtBand.setDescription(srcBand.getDescription());
                sourceProduct.addBand(virtBand);

                final Band targetBand = ProductUtils.copyBand(virtBand.getName(), sourceProduct, targetProduct);
                targetBand.setName(srcBand.getName());
                targetBand.setSourceImage(virtBand.getSourceImage());
            } else {
                final Band targetBand = ProductUtils.copyBand(srcBand.getName(), sourceProduct, targetProduct);
                if(byPass) {
                    targetBand.setSourceImage(srcBand.getSourceImage());
                }
            }
        }
    }

    public void dispose() {
        if(geometry != null && !geometry.isEmpty() && !byPass) {
            final Band[] sourceBands = sourceProduct.getBands();
            for (Band srcBand : sourceBands) {
                if(srcBand.getName().contains(tmpVirtBandName)) {
                    sourceProduct.removeBand(srcBand);
                }
            }
        }
    }

    /**
     * Called by the framework in order to compute the stack of tiles for the given target bands.
     * <p>The default implementation throws a runtime exception with the message "not implemented".</p>
     *
     * @param targetTiles     The current tiles to be computed for each target band.
     * @param targetRectangle The area in pixel coordinates to be computed (same for all rasters in <code>targetRasters</code>).
     * @param pm              A progress monitor which should be used to determine computation cancelation requests.
     * @throws OperatorException if an error occurs during computation of the target rasters.
     */
    @Override
    public void computeTileStack(Map<Band, Tile> targetTiles, Rectangle targetRectangle, ProgressMonitor pm) throws OperatorException {

        try {
            if(dem == null) {
                createDEM();
            }

            final GeoPos geoPos = new GeoPos();
            final PixelPos pixelPos = new PixelPos();
            final GeoCoding targetGeoCoding = targetProduct.getGeoCoding();

            final TileData[] trgTiles = getTargetTiles(targetTiles, targetRectangle, sourceProduct);

            final int maxX = targetRectangle.x + targetRectangle.width - 1;
            final int maxY = targetRectangle.y + targetRectangle.height - 1;
            boolean valid;

            final Tile firstTgtTile = trgTiles[0].targetTile;
            final Tile firstSrcTile = trgTiles[0].srcTile;

            final float demNoDataValue = dem.getDescriptor().getNoDataValue();

            for (int y = targetRectangle.y; y <= maxY; ++y) {
                pixelPos.y = y + 0.5f;
                for (int x = targetRectangle.x; x <= maxX; ++x) {
                    pixelPos.x = x + 0.5f;                    
                    targetGeoCoding.getGeoPos(pixelPos, geoPos);

                    final int tgtIndex = firstTgtTile.getDataBufferIndex(x, y);
                    
                    if(landMask) {
                        if(useSRTM)
                            valid = dem.getElevation(geoPos) == demNoDataValue;
                        else
                            valid = dem.getElevation(geoPos) < seaThreshold;
                    } else {
                        if(useSRTM)
                            valid = dem.getElevation(geoPos) != demNoDataValue;
                        else
                            valid = dem.getElevation(geoPos) > landThreshold;
                    }

                    if(valid) {
                        final int srcIndex = firstSrcTile.getDataBufferIndex(x, y);
                        for(TileData tileData : trgTiles) {
                            tileData.tileDataBuffer.setElemDoubleAt(tgtIndex,
                                    tileData.srcDataBuffer.getElemDoubleAt(srcIndex));
                        }
                    } else {
                        for(TileData tileData : trgTiles) {
                            tileData.tileDataBuffer.setElemDoubleAt(tgtIndex, tileData.noDataValue);
                        }
                    }
                }
            }

        } catch (Throwable e) {
            OperatorUtils.catchOperatorException(getId(), e);
        } finally {
            pm.done();
        }
    }

    private synchronized void createDEM() {
        if(dem != null) return;

        final ElevationModelRegistry elevationModelRegistry = ElevationModelRegistry.getInstance();
        ElevationModelDescriptor demDescriptor = elevationModelRegistry.getDescriptor("ACE2_5Min");
        if(useSRTM) {
            demDescriptor = elevationModelRegistry.getDescriptor("SRTM 3Sec GeoTiff");
        }
        if (demDescriptor.isInstallingDem()) {
              throw new OperatorException("The DEM is currently being installed.");
        }
        dem = demDescriptor.createDem(ResamplingFactory.createResampling(ResamplingFactory.NEAREST_NEIGHBOUR_NAME));
    }

    private TileData[] getTargetTiles(final Map<Band, Tile> targetTiles, final Rectangle targetRectangle,
                                             final Product srcProduct) {
        final ArrayList<TileData> trgTileList = new ArrayList<TileData>();
        final Set<Band> keySet = targetTiles.keySet();
        for(Band targetBand : keySet) {

            final TileData td = new TileData();
            td.targetTile = targetTiles.get(targetBand);
            td.srcTile = getSourceTile(srcProduct.getBand(targetBand.getName()), targetRectangle);
            td.tileDataBuffer = td.targetTile.getDataBuffer();
            td.srcDataBuffer = td.srcTile.getDataBuffer();
            td.noDataValue = targetBand.getNoDataValue();
            trgTileList.add(td);
        }
        return trgTileList.toArray(new TileData[trgTileList.size()]);
    }

    private static class TileData {
        Tile targetTile = null;
        Tile srcTile = null;
        ProductData tileDataBuffer = null;
        ProductData srcDataBuffer = null;
        double noDataValue = 0;
    }

    /**
     * Operator SPI.
     */
    public static class Spi extends OperatorSpi {

        public Spi() {
            super(CreateLandMaskOp.class);
            setOperatorUI(CreateLandMaskOpUI.class);
        }
    }
}