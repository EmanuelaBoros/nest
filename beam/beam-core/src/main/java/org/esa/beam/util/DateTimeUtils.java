/*
 * $Id: DateTimeUtils.java,v 1.2 2009-07-07 00:27:41 lveci Exp $
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
package org.esa.beam.util;

import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.Date;


/**
 * This utility class provides some date/time related methods.
 *
 * @author Norman Fomferra
 * @version $Revision: 1.2 $ $Date: 2009-07-07 00:27:41 $
 * @see java.util.Date
 */
public class DateTimeUtils {

    /**
     * An ISO 8601 date/time format.
     */
    public static final SimpleDateFormat ISO_8601_FORMAT = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss.S");

    /**
     * The number of days from Jan 1, 4713 BC (Proleptic Julian) to 1/1/1970 AD (Gregorian).  1/1/1970 is time zero for
     * a java.util.Date.
     */
    public static final double JD_OFFSET = 2440588.0;
    /**
     * The Modified Julian Day (MJD) gives the number of days since midnight on November 17, 1858. This date corresponds
     * to <code>MJD_OFFSET = 2400000.5</code> days after day zero of the Julian calendar.
     */
    public static final double MJD_OFFSET = 2400000.5;


    /**
     * The number of hours per day.
     */
    public static final double HOURS_PER_DAY = 24.0;
    /**
     * The number of seconds per day.
     */
    public static final double SECONDS_PER_DAY = 3600.0 * HOURS_PER_DAY;
    /**
     * The number of milli-seconds per day.
     */
    public static final double MILLIS_PER_DAY = 1000.0 * SECONDS_PER_DAY;
    /**
     * The number of micro-seconds per day.
     */
    public static final double MICROS_PER_DAY = 1000.0 * MILLIS_PER_DAY;


    /**
     * Converts a julian day (JD) to a modified julian day (MJD) value.
     *
     * @param jd the julian day
     *
     * @return the modified julian day
     */
    public static double jdToMJD(double jd) {
        return jd - MJD_OFFSET;
    }

    /**
     * Converts a modified julian day (MJD) to a julian day (JD) value.
     *
     * @param mjd the modified julian day
     *
     * @return the julian day
     */
    public static double mjdToJD(double mjd) {
        return MJD_OFFSET + mjd;
    }

    /**
     * Converts a julian day (JD) to a UTC date/time value.
     * <p/>
     * <p><i>Important note:</i> Due to the limitations of {@link java.util.Date java.util.Date} this method does not
     * take leap seconds into account.
     *
     * @param jd the julian day
     *
     * @return the UTC date/time
     */
    public static Date jdToUTC(double jd) {
        long millis = Math.round((jd - JD_OFFSET) * MILLIS_PER_DAY);
        return new Date(millis);
    }

    /**
     * Converts a UTC date/time value to a julian day (JD).
     * <p/>
     * <p><i>Important note:</i> Due to the limitations of {@link java.util.Date java.util.Date} this method does not
     * take leap seconds into account.
     *
     * @param utc the UTC date/time, if <code>null</code> the current time is converted
     *
     * @return the julian day
     */
    public static double utcToJD(Date utc) {
        long millis = utc != null ? utc.getTime() : System.currentTimeMillis();
        return JD_OFFSET + millis / MILLIS_PER_DAY;
    }

    /**
     * Converts a UTC date/time value to a string. The method uses the ISO 8601 date/time format <code>YYYY-MM-DD
     * hh:mm:ss.S</code>
     * <p/>
     * <p><i>Important note:</i> Due to the limitations of {@link java.util.Date java.util.Date} this method does not
     * take leap seconds into account.
     *
     * @param utc the UTC date/time value
     *
     * @return the UTC date/time string
     */
    public static String utcToString(Date utc) {
        return ISO_8601_FORMAT.format(utc != null ? utc : new Date());
    }

    /**
     * Converts a UTC date/time string to a UTC date/time value. The method uses the ISO 8601 date/time format
     * <code>YYYY-MM-DD hh:mm:ss.S</code>
     * <p/>
     * <p><i>Important note:</i> Due to the limitations of {@link java.util.Date java.util.Date} this method does not
     * take leap seconds into account.
     *
     * @param utc the UTC date/time string
     */
    public static Date stringToUTC(String utc) throws ParseException {
        return ISO_8601_FORMAT.parse(utc);
    }
}
