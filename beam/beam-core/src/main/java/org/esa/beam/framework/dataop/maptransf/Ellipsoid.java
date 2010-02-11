/*
 * $Id: Ellipsoid.java,v 1.3 2010-02-11 17:02:24 lveci Exp $
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
package org.esa.beam.framework.dataop.maptransf;

/**
 * Represents an ellipsoid used to approximate the earth's surface.
 * 
 * @deprecated since BEAM 4.7, use {@link org.opengis.referencing.datum.Ellipsoid} instead.
 */
@Deprecated
public class Ellipsoid {

    /**
     * The standard WGS-72 ellipsoid.
     */
    public static final Ellipsoid WGS_72 = new Ellipsoid("WGS-72", 6356750.5, 6378135.0);

    /**
     * The standard WGS-84 ellipsoid.
     */
    public static final Ellipsoid WGS_84 = new Ellipsoid("WGS-84", 6356752.3, 6378137.0);

    /**
     * Bessel ellipsoid (used for Gauss-Krueger Projection)
     */
    public static final Ellipsoid BESSEL = new Ellipsoid("Bessel 1841", 6356079.0, 6377397.2);

    /**
     * The GRS-80 ellipsoid
     */
    public static final Ellipsoid GRS_80 = new Ellipsoid("GRS-80", 6356752.3141, 6378137.0);

    private final String name;
    private final double semiMinor;
    private final double semiMajor;

    public Ellipsoid(String name, double semiMinor, double semiMajor) {
        this.name = name;
        this.semiMinor = semiMinor;
        this.semiMajor = semiMajor;
    }

    public String getName() {
        return name;
    }

    public double getSemiMinor() {
        return semiMinor;
    }

    public double getSemiMajor() {
        return semiMajor;
    }

    @Override
    public Object clone() {
        try {
            return super.clone();
        } catch (CloneNotSupportedException e) {
            throw new IllegalStateException(e);
        }
    }

    @Override
    public boolean equals(Object obj) {
        if (obj instanceof Ellipsoid) {
            final Ellipsoid that = (Ellipsoid) obj;
            return this.name.equals(that.name) &&
                   this.semiMajor == that.semiMajor &&
                   this.semiMinor == that.semiMinor;
        }

        return false;
    }

    @Override
    public int hashCode() {
        int result = 17;

        result = 31 * result + name.hashCode();
        result = 31 * result + hashCodeDouble(semiMajor);
        result = 31 * result + hashCodeDouble(semiMinor);

        return result;
    }

    private static int hashCodeDouble(double d) {
        final long l = Double.doubleToLongBits(d);
        return (int) (l ^ (l >>> 32));
    }
}
