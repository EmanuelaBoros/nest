/*
 * $Id: IntegerizedSinusoidalDescriptor.java,v 1.2 2010-02-11 17:02:24 lveci Exp $
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
import org.esa.beam.util.Guardian;

/**
 * <p><i>Note that this class is not yet public API and may change in future releases.</i></p>
 *
 * @author Norman Fomferra
 * @version $Revision: 1.2 $ $Date: 2010-02-11 17:02:24 $
 * 
 * @deprecated since BEAM 4.7, use geotools {@link org.geotools.referencing.operation.MathTransformProvider} instead.
 */
@Deprecated
public class IntegerizedSinusoidalDescriptor implements MapTransformDescriptor {

    public static final String TYPE_ID = "ISEAG";
    public static final String NAME = "Integerized Sinusoidal Equal Area Grid";
    public static final String MAP_UNIT = "units";
    public static final int DEFAULT_ROW_COUNT_VALUE = 2160;
    private static final Parameter[] PARAMETERS = new Parameter[] {new Parameter("rowCount", DEFAULT_ROW_COUNT_VALUE)};
    private static final double[] PARAMETER_DEFAULT_VALUES = new double[] {DEFAULT_ROW_COUNT_VALUE};

    public void registerProjections() {
        MapProjectionRegistry.registerProjection(new MapProjection(getName(), createTransform(null), false));
    }

    public IntegerizedSinusoidalDescriptor() {
    }

    public String getTypeID() {
        return TYPE_ID;
    }

    /**
     * Gets a descriptive name for this map transformation descriptor, e.g. "Transverse Mercator".
     */
    public String getName() {
        return NAME;
    }

    public String getMapUnit() {
        return MAP_UNIT;
    }

    /**
     * Gets the default parameter values for this map transform.
     */
    public double[] getParameterDefaultValues() {
        final double[] values = new double[PARAMETER_DEFAULT_VALUES.length];
    	System.arraycopy(PARAMETER_DEFAULT_VALUES, 0, values, 0, values.length);

        return values;
    }

    public Parameter[] getParameters() {
        return PARAMETERS;
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
        return new DefaultMapTransformUI(transform);
    }

    public MapTransform createTransform(double[] parameterValues) {
        if (parameterValues == null) {
            parameterValues = getParameterDefaultValues();
        }
        return new MT(parameterValues);
    }

    private class MT implements MapTransform {

        private final ISEAG.RC _rcTemp = new ISEAG.RC();
        private final ISEAG.LL _llTemp = new ISEAG.LL();
        private ISEAG _grid;

        private MT(double[] parameterValues) {
            this((int) parameterValues[0]);
        }

        private MT(int rowCount) {
            Guardian.assertGreaterThan("rowCount", rowCount, 3);
            _grid = new ISEAG(rowCount);
        }

        public MapTransformDescriptor getDescriptor() {
            return IntegerizedSinusoidalDescriptor.this;
        }

        public double[] getParameterValues() {
            return new double[]{_grid.getRowCount()};
        }

        /**
         * Forward project geographical co-ordinates into map co-ordinates.
         */
        public Point2D forward(GeoPos geoPoint, Point2D mapPoint) {
            if (mapPoint == null) {
                mapPoint = new Point2D.Float();
            }
            _grid.ll2rc(geoPoint.lat, geoPoint.lon, _rcTemp);
            final int ncols = _grid.getColumnCount(_rcTemp.row);
            final double x = _rcTemp.col - 0.5 * ncols;
            final double y = _rcTemp.row;
            mapPoint.setLocation(x,  y);
            return mapPoint;
        }

        /**
         * Inverse project map co-ordinates into geographical co-ordinates.
         */
        public GeoPos inverse(Point2D mapPoint, GeoPos geoPoint) {
            if (geoPoint == null) {
                geoPoint = new GeoPos();
            }
            final int ncols = _grid.getColumnCount(_rcTemp.row);
            final double x = mapPoint.getX();
            final double y = mapPoint.getY();
            final int col = (int)Math.floor(x + 0.5 * ncols);
            final int row = (int)Math.floor(y);
            _grid.rc2ll(row, col, _llTemp);
            geoPoint.lon = (float) _llTemp.lon;
            geoPoint.lat = (float) _llTemp.lat;
            return geoPoint;
        }

        public MapTransform createDeepClone() {
            return new MT(_grid.getRowCount());
        }
    }
}
