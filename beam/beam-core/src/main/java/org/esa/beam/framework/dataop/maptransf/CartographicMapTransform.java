/*
 * $Id: CartographicMapTransform.java,v 1.3 2010-02-11 17:02:24 lveci Exp $
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

import org.esa.beam.framework.datamodel.GeoPos;

import java.awt.geom.Point2D;

/**
 * An abstract base class for cartographic map-transformations.
 * 
 * @deprecated since BEAM 4.7, use geotools {@link org.geotools.referencing.operation.projection.MapProjection} instead.
 */
@Deprecated
public abstract class CartographicMapTransform implements MapTransform {

    /**
     * The central meridian value in degree.
     */
    protected final float _centralMeridian;
    /**
     * Semi-major parameter of ellipsoid (map scaling factor) in meter.
     */
    protected final double _a;
    /**
     * Inverse semi-major parameter.
     */
    protected final double _invA;
    /**
     * The false easting (map x-offset) in map units.
     */
    protected final double _x0;
    /**
     * The false northing (map x-offset) in map units.
     */
    protected final double _y0;

    /**
     * Constructs a new instance.
     *
     * @param centralMeridian central meridian in degree
     * @param falseEasting    false easting (map x-offset) in map units
     * @param falseNorthing   false northing (map y-offset) in map units
     * @param semiMajor       semi-major parameter of ellipsoid (map scaling factor) in map units
     */
    protected CartographicMapTransform(double centralMeridian,
                                       double falseEasting,
                                       double falseNorthing,
                                       double semiMajor) {
        _centralMeridian = (float) centralMeridian;
        _x0 = falseEasting;
        _y0 = falseNorthing;
        _a = semiMajor;
        _invA = 1.0 / semiMajor;
    }

    /**
     * Gets the central meridian parameter.
     *
     * @return the central meridian
     */
    public double getCentralMeridian() {
        return _centralMeridian;
    }

    /**
     * Gets the map scaling factor parameter.
     *
     * @return the map scaling factor
     */
    public double getSemiMajor() {
        return _a;
    }

    /**
     * Gets the map scaling factor parameter.
     *
     * @return the map scaling factor
     */
    public double getInverseSemiMajor() {
        return _invA;
    }

    /**
     * Gets the map X-offset (false easting) parameter.
     *
     * @return the map X-offset
     */
    public double getFalseEasting() {
        return _x0;
    }

    /**
     * Gets the map Y-offset (false northing) parameter.
     *
     * @return the map Y-offset
     */
    public double getFalseNorthing() {
        return _y0;
    }

    /**
     * Forward project geographical co-ordinates into map co-ordinates.
     *
     * @param geoPoint the source position in lat/lon
     * @param mapPoint the target map position in x/y (,ight be null, then a new object is created
     */
    @Override
    public Point2D forward(GeoPos geoPoint, Point2D mapPoint) {
        if (mapPoint == null) {
            mapPoint = new Point2D.Double();
        }

        float lon0 = geoPoint.getLon() - _centralMeridian;
        if (lon0 < -180.0) {
            lon0 = (float) (360.0 + lon0);
        }
        if (lon0 > 180.0) {
            lon0 = (float) (lon0 - 360.0);
        }

        mapPoint = forward_impl(geoPoint.getLat(),
                                lon0,
                                mapPoint);
        mapPoint.setLocation(_a * mapPoint.getX() + _x0,
                             _a * mapPoint.getY() + _y0);
        return mapPoint;
    }

    /**
     * Inverse project map co-ordinates into geographical co-ordinates.
     *
     * @param mapPoint the source location in x/y map coordinates
     * @param geoPoint the target position in lat/lon (might be null, then a new object is created)
     */
    @Override
    public GeoPos inverse(Point2D mapPoint, GeoPos geoPoint) {
        if (geoPoint == null) {
            geoPoint = new GeoPos();
        }
        geoPoint = inverse_impl((mapPoint.getX() - _x0) * _invA,
                                (mapPoint.getY() - _y0) * _invA,
                                geoPoint);
        geoPoint.setLocation(geoPoint.getLat(),
                             geoPoint.getLon() + _centralMeridian);

        return geoPoint;
    }

    ///////////////////////////////////////////////////////////////////////////
    /////// END OF PUBLIC
    ///////////////////////////////////////////////////////////////////////////

    /**
     * Worker method to be overridden by derived class. Performs the pure transformation. Prescaling, northing, easting
     * etc is calculated in this class.
     *
     * @param lat      latitude of source location
     * @param lon      longitude of source location
     * @param mapPoint point on the map
     *
     * @return the map coordinate
     */
    abstract protected Point2D forward_impl(float lat, float lon, Point2D mapPoint);

    /**
     * Worker method to be overridden by derived class. Performs the pure transformation. Prescaling, northing, easting
     * etc is calculated in this class.
     * <p/>
     * <p>Should be overridden in order to delegate to <code>{@link #inverse_impl(double, double,
     * org.esa.beam.framework.datamodel.GeoPos)}</code> if transformation is performed is in 64-bit accuracy. Override
     * <code>{@link #inverse_impl(double, double, org.esa.beam.framework.datamodel.GeoPos)}</code> instead in order to
     * perform the actual transformation.
     *
     * @param geoPoint point on the earth's surface
     * @param x        map x coordinate
     * @param y        map y coordinate
     *
     * @return the geodetic co-ordinate
     */
    abstract protected GeoPos inverse_impl(float x, float y, GeoPos geoPoint);

    /**
     * Worker method to be overridden by derived class. Performs the pure transformation. Prescaling, northing, easting
     * etc is calculated in this class.
     * <p/>
     * <p>Should be overridden in order to perform transformation in 64-bit accuracy.
     * <p/>
     * <p>The default implementation simple returns <code>inverse_impl((float)x, (float)y, geoPoint)</code>.
     *
     * @param geoPoint point on the earth's surface
     * @param x        map x coordinate
     * @param y        map y coordinate
     *
     * @return the geodetic co-ordinate
     */
    protected GeoPos inverse_impl(double x, double y, GeoPos geoPoint) {
        return inverse_impl((float) x, (float) y, geoPoint);
    }
}
