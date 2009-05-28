/*
 * $Id: Datum.java,v 1.2 2009-05-27 21:09:23 lveci Exp $
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
 * Represents a geodetic datum. Geodetic datums define the size and shape of the earth and the origin and orientation of
 * the coordinate systems used to map the earth.
 */
public class Datum implements Cloneable {

    /**
     * The standard WGS-72 datum.
     */
    public static final Datum WGS_72 = new Datum("WGS-72", Ellipsoid.WGS_72, 0.0, 0.0, 5.0);

    /**
     * The standard WGS-84 datum.
     */
    public static final Datum WGS_84 = new Datum("WGS-84", Ellipsoid.WGS_84, 0.0, 0.0, 0.0);

    /**
     * The ITRF-97 datum.
     */
    public static final Datum ITRF_97 = new Datum("ITRF-97", Ellipsoid.GRS_80, 0.0, 0.0, 0.0);

    private final String name;
    private final Ellipsoid ellipsoid;
    private final double dx;
    private final double dy;
    private final double dz;

    public Datum(String name, Ellipsoid ellipsoid, double dx, double dy, double dz) {
        this.name = name;
        this.ellipsoid = ellipsoid;
        this.dx = dx;
        this.dy = dy;
        this.dz = dz;
    }

    public String getName() {
        return name;
    }

    public Ellipsoid getEllipsoid() {
        return ellipsoid;
    }


    public double getDX() {
        return dx;
    }

    public double getDY() {
        return dy;
    }

    public double getDZ() {
        return dz;
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
        if (obj instanceof Datum) {
            final Datum that = (Datum) obj;
            return this.name.equals(that.name) &&
                   this.dx == that.dx &&
                   this.dy == that.dy &&
                   this.dz == that.dz &&
                   this.ellipsoid.equals(that.ellipsoid);
        }

        return false;
    }

    @Override
    public int hashCode() {
        int result = 17;

        result = 31 * result + name.hashCode();
        result = 31 * result + hashCodeDouble(dx);
        result = 31 * result + hashCodeDouble(dy);
        result = 31 * result + hashCodeDouble(dz);
        result = 31 * result + ellipsoid.hashCode();

        return result;
    }

    private static int hashCodeDouble(double d) {
        final long l = Double.doubleToLongBits(d);
        return (int) (l ^ (l >>> 32));
    }
}
