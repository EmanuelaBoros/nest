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
package org.esa.nest.dataio.dem.srtm3_geotiff;

import org.esa.beam.framework.dataop.dem.AbstractElevationModelDescriptor;
import org.esa.beam.framework.dataop.dem.ElevationModel;
import org.esa.beam.framework.dataop.maptransf.Datum;
import org.esa.beam.framework.dataop.resamp.Resampling;
import org.esa.beam.util.SystemUtils;
import org.esa.nest.util.Settings;

import java.io.File;
import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URL;

public class SRTM3GeoTiffElevationModelDescriptor extends AbstractElevationModelDescriptor {

    public static final String NAME = "SRTM 3Sec GeoTiff";
    public static final String DB_FILE_SUFFIX = ".TIF";
    public static final String ARCHIVE_URL_PATH = SystemUtils.BEAM_HOME_PAGE + "data/ACE.zip";
    public static final int NUM_X_TILES = 72;
    public static final int NUM_Y_TILES = 24;
    public static final int DEGREE_RES = 5;
    public static final int PIXEL_RES = 6000;
    public static final int NO_DATA_VALUE = -32768;
    public static final int RASTER_WIDTH = NUM_X_TILES * PIXEL_RES;
    public static final int RASTER_HEIGHT = NUM_Y_TILES * PIXEL_RES;
    public static final Datum DATUM = Datum.WGS_84;

    private File demInstallDir = null;

    public SRTM3GeoTiffElevationModelDescriptor() {
    }

    public String getName() {
        return NAME;
    }

    public Datum getDatum() {
        return DATUM;
    }

    public float getNoDataValue() {
        return NO_DATA_VALUE;
    }

    @Override
    public File getDemInstallDir() {
        if(demInstallDir == null) {
            final String path = Settings.instance().get("DEM/srtm3GeoTiffDEMDataPath");
            demInstallDir = new File(path);
            if(!demInstallDir.exists())
                demInstallDir.mkdirs();
        }
        return demInstallDir;
    }

    public boolean isDemInstalled() {
        return true;
    }

    public URL getDemArchiveUrl() {
        try {
            return new URL(ARCHIVE_URL_PATH);
        } catch (MalformedURLException e) {
            throw new IllegalStateException("MalformedURLException not expected: " + ARCHIVE_URL_PATH);
        }
    }

    @Deprecated
    public ElevationModel createDem() {
        try {
            return new SRTM3GeoTiffElevationModel(this, Resampling.BILINEAR_INTERPOLATION);
        } catch (Exception e) {
            return null;
        }
    }

    public ElevationModel createDem(Resampling resamplingMethod) {
        try {
            return new SRTM3GeoTiffElevationModel(this, resamplingMethod);
        } catch (Exception e) {
            return null;
        }
    }

    public static String createTileFilename(int tileX, int tileY) {
        String name = "srtm_";
        if(tileX < 10)
            name += "0" + tileX;
        else
            name += tileX;
        name += '_';
        if(tileY < 10)
            name += "0" + tileY;
        else
            name += tileY;
        return name + ".zip";
    }

}