/*
 * $Id: CeosHelper.java,v 1.1 2008-01-04 16:23:10 lveci Exp $
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
package org.esa.nest.dataio.ceos;

import org.esa.nest.dataio.ceos.records.FilePointerRecord;
import org.esa.nest.dataio.ceos.records.TextRecord;
import org.esa.nest.dataio.ceos.records.VolumeDescriptorRecord;
import org.esa.beam.framework.datamodel.ProductData;

import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.Calendar;

public class CeosHelper {

    private static final String VOLUME_FILE_PREFIX = "VDF_";
    private static final String LEADER_FILE_PREFIX = "LEA_";
    private static final String IMAGE_FILE_PREFIX = "DAT_";
    private static final String TRAILER_FILE_PREFIX = "NUL_";

    public static File getVolumeFile(final File baseDir) throws IOException {
        final File[] files = baseDir.listFiles(new FilenameFilter() {
            public boolean accept(final File dir, final String name) {
                return name.startsWith(VOLUME_FILE_PREFIX);
            }
        });
        if (files == null || files.length < 1) {
            throw new IOException("No volume descriptor file found in directory:\n"
                                  + baseDir.getPath());
        }
        if (files.length > 1) {
            throw new IOException("Multiple volume descriptor files found in directory:\n"
                                  + baseDir.getPath());
        }
        return files[0];
    }

    public static FilePointerRecord[] readFilePointers(final VolumeDescriptorRecord vdr) throws
                                                                                         IllegalCeosFormatException,
                                                                                         IOException {
        final int numFilePointers = vdr.getNumberOfFilepointerRecords();
        final CeosFileReader reader = vdr.getReader();
        reader.seek(vdr.getRecordLength());
        final FilePointerRecord[] filePointers = new FilePointerRecord[numFilePointers];
        for (int i = 0; i < numFilePointers; i++) {
            filePointers[i] = new FilePointerRecord(reader);
        }
        return filePointers;
    }

    public static String getLeaderFileName(final TextRecord textRecord) {
        return LEADER_FILE_PREFIX + getProductName(textRecord);
    }

    public static String getTrailerFileName(final TextRecord textRecord) {
        return TRAILER_FILE_PREFIX + getProductName(textRecord);
    }

    public static String getImageFileName(final TextRecord textRecord, final String ccd) {
        if (ccd != null && ccd.trim().length() > 0) {
            return IMAGE_FILE_PREFIX + "0" + ccd + "-" + getProductName(textRecord);
        } else {
            return IMAGE_FILE_PREFIX + getProductName(textRecord);
        }
    }

    public static String getProductName(final TextRecord textRecord) {
        return textRecord.getSceneID() + "-" + textRecord.getProductID();
    }

    public static ProductData.UTC createUTCDate(final int year, final int dayOfYear, final int millisInDay) {
        final Calendar calendar = ProductData.UTC.createCalendar();

        calendar.set(Calendar.YEAR, year);
        calendar.set(Calendar.DAY_OF_YEAR, dayOfYear);
        calendar.add(Calendar.MILLISECOND, millisInDay);

        return ProductData.UTC.create(calendar.getTime(), 0);
    }

    public static double[] sortToFXYSumOrder(final double[] coeffs) {
        final double[] newOrder = new double[coeffs.length];
        newOrder[0] = coeffs[0];
        newOrder[1] = coeffs[1];
        newOrder[2] = coeffs[2];
        newOrder[3] = coeffs[4];
        newOrder[4] = coeffs[3];
        newOrder[5] = coeffs[5];
        newOrder[6] = coeffs[8];
        newOrder[7] = coeffs[6];
        newOrder[8] = coeffs[7];
        newOrder[9] = coeffs[9];

        return newOrder;
    }

    public static double[] convertLongToDouble(final long[] longs) {
        final double[] doubles = new double[longs.length];
        for (int i = 0; i < longs.length; i++) {
            doubles[i] = Double.longBitsToDouble(longs[i]);
        }
        return doubles;
    }
}
