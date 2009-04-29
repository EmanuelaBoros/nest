/*
 * $Id: GETASSE30ElevationModelDescriptor.java,v 1.3 2009-04-29 15:00:08 lveci Exp $
 *
 * Copyright (C) 2002 by Brockmann Consult (info@brockmann-consult.de)
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
package org.esa.beam.dataio.getasse30;

import org.esa.beam.framework.dataop.dem.AbstractElevationModelDescriptor;
import org.esa.beam.framework.dataop.dem.ElevationModel;
import org.esa.beam.framework.dataop.maptransf.Datum;
import org.esa.beam.framework.dataop.resamp.Resampling;
import org.esa.beam.util.SystemUtils;

import java.io.File;
import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URL;

public class GETASSE30ElevationModelDescriptor extends AbstractElevationModelDescriptor {

    public static final String NAME = "GETASSE30";
    public static final String DB_FILE_SUFFIX = ".GETASSE30";
    public static final String ARCHIVE_URL_PATH = SystemUtils.BEAM_HOME_PAGE + "data/GETASSE30.zip";
    public static final int NUM_X_TILES = 24;
    public static final int NUM_Y_TILES = 12;
    public static final int DEGREE_RES = 15;
    public static final int PIXEL_RES = 1800;
    public static final int NO_DATA_VALUE = -9999;
    public static final int RASTER_WIDTH = NUM_X_TILES * PIXEL_RES;
    public static final int RASTER_HEIGHT = NUM_Y_TILES * PIXEL_RES;
    public static final Datum DATUM = Datum.WGS_84;

    public GETASSE30ElevationModelDescriptor() {
    }

    @Override
    public String getName() {
        return NAME;
    }

    @Override
    public Datum getDatum() {
        return DATUM;
    }

    @Override
    public float getNoDataValue() {
        return NO_DATA_VALUE;
    }

    @Override
    public boolean isDemInstalled() {
        final File file = getTileFile(-180, -90);   // todo (nf) - check all tiles
        return file.canRead();
    }

    @Override
    public URL getDemArchiveUrl() {
        try {
            return new URL(ARCHIVE_URL_PATH);
        } catch (MalformedURLException e) {
            throw new IllegalStateException("MalformedURLException not expected: " + ARCHIVE_URL_PATH);
        }
    }

    @Override
    @Deprecated
    public ElevationModel createDem() {
        return createDem(Resampling.BILINEAR_INTERPOLATION);
    }

    @Override
    public ElevationModel createDem(Resampling resampling) {
        try {
            return new GETASSE30ElevationModel(this, resampling);
        } catch (IOException e) {
            return null;
        }
    }

    public File getTileFile(int minLon, int minLat) {
        return new File(getDemInstallDir(), createTileFilename(minLat, minLon));
    }

    public String createTileFilename(int minLat, int minLon) {
        String latString = minLat < 0 ? Math.abs(minLat) + "S" : minLat + "N";
        while (latString.length() < 3) {
            latString = "0" + latString;
        }
        String lonString = minLon < 0 ? Math.abs(minLon) + "W" : minLon + "E";
        while (lonString.length() < 4) {
            lonString = "0" + lonString;
        }
        return latString + lonString + DB_FILE_SUFFIX;
    }

}
