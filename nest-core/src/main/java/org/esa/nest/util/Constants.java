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
package org.esa.nest.util;

public final class Constants
{
    public static final double lightSpeed = 299792458.0; //  m / s
    public static final double halfLightSpeed = lightSpeed / 2.0;
    public static final double lightSpeedInMetersPerDay = Constants.lightSpeed * 86400.0;

    public static final double semiMajorAxis = 6378137.0;      // in m, WGS84 semi-major axis of Earth
    public static final double semiMinorAxis = 6356752.314245; // in m, WGS84 semi-minor axis of Earth

    public static final double MeanEarthRadius = 6371008.7714; // in m (WGS84)

    public static final double oneMillion = 1000000.0;
    public static final double oneBillion = 1000000000.0;

    public static final double TWO_PI = 2.0*Math.PI;

    public static final double EPS = 1e-15;

    private Constants()
    {
    }
}