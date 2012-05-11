/*
 * Copyright (C) 2011 by Array Systems Computing Inc. http://www.array.ca
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
package org.esa.nest.dataio.dem.srtm3_geotiff;

import org.esa.beam.framework.dataio.ProductReaderPlugIn;
import org.esa.beam.framework.datamodel.GeoPos;
import org.esa.beam.framework.datamodel.PixelPos;
import org.esa.beam.framework.dataop.resamp.Resampling;
import org.esa.nest.dataio.dem.BaseElevationModel;
import org.esa.nest.dataio.dem.ElevationFile;

import java.io.File;

public final class SRTM3GeoTiffElevationModel extends BaseElevationModel {

    private static final ProductReaderPlugIn productReaderPlugIn = getReaderPlugIn("GeoTIFF");

    public SRTM3GeoTiffElevationModel(final SRTM3GeoTiffElevationModelDescriptor descriptor, final Resampling resamplingMethod) {
        super(descriptor, resamplingMethod);

        setMaxCacheSize(12);
    }

    @Override
    public PixelPos getIndex(final GeoPos geoPos) {
        final float pixelY = RASTER_HEIGHT - (geoPos.lat + 60.0f) / DEGREE_RES_BY_NUM_PIXELS_PER_TILE;
        final float pixelX = (geoPos.lon + 180.0f) / DEGREE_RES_BY_NUM_PIXELS_PER_TILE; 
        return new PixelPos(pixelX, pixelY);
    }

    @Override
    public GeoPos getGeoPos(final PixelPos pixelPos) {
        final float pixelLat = (RASTER_HEIGHT - pixelPos.y) * DEGREE_RES_BY_NUM_PIXELS_PER_TILE - 60.0f;
        final float pixelLon = pixelPos.x * DEGREE_RES_BY_NUM_PIXELS_PER_TILE - 180.0f;
        return new GeoPos(pixelLat, pixelLon);
    }

    @Override
    protected void createElevationFile(final ElevationFile[][] elevationFiles,
                                                final int x, final int y, final File demInstallDir) {
        final String fileName = descriptor.createTileFilename(x + 1, y + 1);
        final File localFile = new File(demInstallDir, fileName);
        elevationFiles[x][y] = new SRTM3GeoTiffFile(this, localFile, productReaderPlugIn.createReaderInstance());
    }
}