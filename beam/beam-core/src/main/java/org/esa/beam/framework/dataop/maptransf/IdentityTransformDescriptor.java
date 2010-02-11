/*
 * $Id: IdentityTransformDescriptor.java,v 1.2 2010-02-11 17:02:24 lveci Exp $
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

import java.awt.geom.Point2D;

import org.esa.beam.framework.datamodel.GeoPos;
import org.esa.beam.framework.param.Parameter;

/**
 * The descriptor for a map transformation which represents the identity transformation.
 *
 * @author Norman Fomferra
 * @version $Revision: 1.2 $ $Date: 2010-02-11 17:02:24 $
 * 
 * @deprecated since BEAM 4.7, use geotools {@link org.geotools.referencing.operation.projection.PlateCarree.Provider} instead.
 */
@Deprecated
public class IdentityTransformDescriptor implements MapTransformDescriptor {

    public static final String TYPE_ID = "Identity";
    public static final String NAME = "Geographic Lat/Lon";
    public static final String MAP_UNIT = "degree";
    public static final Parameter[] PARAMETERS = new Parameter[] {};
    public static final double[] PARAMETER_DEFAULT_VALUES = new double[] {};

    public IdentityTransformDescriptor() {
    }

    public void registerProjections() {
        MapProjectionRegistry.registerProjection(new MapProjection(getName(), createTransform(null), true));
    }

    /**
     * Gets a descriptive name for this map transformation descriptor.
     */
    public String getName() {
        return NAME;
    }

    public String getTypeID() {
        return TYPE_ID;
    }

    public String getMapUnit() {
        return MAP_UNIT;
    }

    public Parameter[] getParameters() {
        return new Parameter[] {};
    }

    /**
     * Gets the default parameter values for this map transform.
     */
    public double[] getParameterDefaultValues() {
        return new double[] {};
    }

    /**
     * Tests if a user interface is available. Returns <code>false</code> because a user interface is not available for
     * this descriptor.
     *
     * @return always <code>false</code>
     */
    public boolean hasTransformUI() {
        return false;
    }

    /**
     * Gets a user interface for editing the transformation properties of a map projection. Returns <code>null</code>
     * because a user interface is not available for this descriptor.
     *
     * @param transform ignored
     *
     * @return always <code>null</code>
     */
    public MapTransformUI getTransformUI(MapTransform transform) {
        return null;
    }

    public MapTransform createTransform(double[] parameterValues) {
        return new IMT();
    }

    private class IMT implements MapTransform {

        private IMT() {
        }

        public MapTransformDescriptor getDescriptor() {
            return IdentityTransformDescriptor.this;
        }

        public double[] getParameterValues() {
            return new double[] {};
        }

        /**
         * Forward project geographical co-ordinates into map co-ordinates.
         */
        public Point2D forward(GeoPos geoPoint, Point2D mapPoint) {
            if (mapPoint != null) {
                mapPoint.setLocation(geoPoint.lon, geoPoint.lat);
            } else {
                mapPoint = new Point2D.Double(geoPoint.lon, geoPoint.lat);
            }
            return mapPoint;
        }

        /**
         * Inverse project map co-ordinates into geographical co-ordinates.
         */
        public GeoPos inverse(Point2D mapPoint, GeoPos geoPoint) {
            if (geoPoint != null) {
                geoPoint.setLocation((float) mapPoint.getY(), (float) mapPoint.getX());
            } else {
                geoPoint = new GeoPos((float) mapPoint.getY(), (float) mapPoint.getX());
            }
            return geoPoint;
        }

        public MapTransform createDeepClone() {
            return new IMT();
        }
    }
}
