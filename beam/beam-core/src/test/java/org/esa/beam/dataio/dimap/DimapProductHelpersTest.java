/*
 * $Id: DimapProductHelpersTest.java,v 1.3 2009-09-01 20:27:12 lveci Exp $
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
package org.esa.beam.dataio.dimap;

import com.bc.ceres.core.ProgressMonitor;
import junit.framework.TestCase;
import org.esa.beam.framework.datamodel.Band;
import org.esa.beam.framework.datamodel.CrsGeoCoding;
import org.esa.beam.framework.datamodel.FXYGeoCoding;
import org.esa.beam.framework.datamodel.GeoCoding;
import org.esa.beam.framework.datamodel.MapGeoCoding;
import org.esa.beam.framework.datamodel.PixelGeoCoding;
import org.esa.beam.framework.datamodel.Product;
import org.esa.beam.framework.datamodel.ProductData;
import org.esa.beam.framework.dataop.maptransf.Datum;
import org.esa.beam.framework.dataop.maptransf.Ellipsoid;
import org.esa.beam.framework.dataop.maptransf.LambertConformalConicDescriptor;
import org.esa.beam.framework.dataop.maptransf.MapInfo;
import org.esa.beam.framework.dataop.maptransf.MapProjection;
import org.esa.beam.framework.dataop.maptransf.MapTransform;
import org.esa.beam.framework.dataop.maptransf.MapTransformDescriptor;
import org.esa.beam.framework.dataop.resamp.Resampling;
import org.esa.beam.util.ArrayUtils;
import org.esa.beam.util.StringUtils;
import org.esa.beam.util.SystemUtils;
import org.esa.beam.util.math.FXYSum;
import org.geotools.referencing.CRS;
import org.jdom.Document;
import org.opengis.referencing.crs.CoordinateReferenceSystem;

import java.awt.geom.AffineTransform;
import java.awt.geom.Rectangle2D;
import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.util.Arrays;
import java.util.Vector;

public class DimapProductHelpersTest extends TestCase {

    private static final String LS = SystemUtils.LS;
    private static final String _projectionName = "ProjectionName";
    private static final String _mapUnit = "mapUnit";
    private static final String _ellipsoidName = "EllipsoidName";
    private static final double _semiMinor = 1234.0;
    private static final double _semiMajor = 5678.0;
    private static final double[] _expValues = new double[]{_semiMajor, _semiMinor, 15, 16, 17, 18, 19};
    private static final String _datumName = "DatumName";
    private static final float _pixelX = 3.2f;
    private static final float _pixelY = 4.3f;
    private static final float _orientation = 11.12f;
    private static final float _easting = 5.4f;
    private static final float _northing = 6.5f;
    private static final float _pixelSizeX = 7.6f;
    private static final float _pixelSizeY = 8.7f;
    private static final double _noDataValue = 1.0e-3;
    private static final boolean _orthorectified = false;
    private static final String _elevModelName = "DEM";
    private static final String _typeId = LambertConformalConicDescriptor.TYPE_ID;
    private static final int _sceneWidth = 200;
    private static final int _sceneHeight = 300;
    private static final boolean _sceneFitted = true;
    private static final String _resamplingName = Resampling.NEAREST_NEIGHBOUR.getName();
    private static final String[] _paramUnit = LambertConformalConicDescriptor.PARAMETER_UNITS;
    private static final String[] _paramNames = LambertConformalConicDescriptor.PARAMETER_NAMES;

    private static final int[] _notMandatoryLines = new int[]{2, 4, 5, 7, 25, 29, 33, 37, 41, 45, 49};

    private final String[] _xmlMapGeocodingStringStyleV1_4_0 = new String[]{
        "<" + DimapProductConstants.TAG_ROOT + ">" + LS,
        /*  1 */ "    <Coordinate_Reference_System>" + LS,
        /*  2 */ "        <GEO_TABLES version=\"1.0\">CUSTOM</GEO_TABLES>" + LS, // notMandatory
        /*  3 */ "        <Horizontal_CS>" + LS,
        /*  4 */ "            <HORIZONTAL_CS_TYPE>PROJECTED</HORIZONTAL_CS_TYPE>" + LS, // notMandatory
        /*  5 */ "            <HORIZONTAL_CS_NAME>" + _projectionName + "</HORIZONTAL_CS_NAME>" + LS, // notMandatory
        /*  6 */ "            <Geographic_CS>" + LS,
        /*  7 */ "                <GEOGRAPHIC_CS_NAME>" + _projectionName + "</GEOGRAPHIC_CS_NAME>" + LS, // notMandatory
        /*  8 */ "                <Horizontal_Datum>" + LS,
        /*  9 */ "                    <HORIZONTAL_DATUM_NAME>" + _datumName + "</HORIZONTAL_DATUM_NAME>" + LS,
        /* 10 */ "                    <Ellipsoid>" + LS,
        /*  1 */ "                        <ELLIPSOID_NAME>" + _ellipsoidName + "</ELLIPSOID_NAME>" + LS,
        /*  2 */ "                        <Ellipsoid_Parameters>" + LS,
        /*  3 */ "                            <ELLIPSOID_MAJ_AXIS unit=\"M\">" + _semiMajor + "</ELLIPSOID_MAJ_AXIS>" + LS,
        /*  4 */ "                            <ELLIPSOID_MIN_AXIS unit=\"M\">" + _semiMinor + "</ELLIPSOID_MIN_AXIS>" + LS,
        /*  5 */ "                        </Ellipsoid_Parameters>" + LS,
        /*  6 */ "                    </Ellipsoid>" + LS,
        /*  7 */ "                </Horizontal_Datum>" + LS,
        /*  8 */ "            </Geographic_CS>" + LS,
        /*  9 */ "            <Projection>" + LS,
        /* 20 */ "                <NAME>" + _projectionName + "</NAME>" + LS,
        /*  1 */ "                <Projection_CT_Method>" + LS,
        /*  2 */ "                    <PROJECTION_CT_NAME>" + _typeId + "</PROJECTION_CT_NAME>" + LS,
        /*  3 */ "                    <Projection_Parameters>" + LS,
        /*  4 */ "                        <Projection_Parameter>" + LS,
        /*  5 */ "                            <PROJECTION_PARAMETER_NAME>" + _paramNames[0] + "</PROJECTION_PARAMETER_NAME>" + LS, // notMandatory
        /*  6 */ "                            <PROJECTION_PARAMETER_VALUE unit=\"" + _paramUnit[0] + "\">" + _expValues[0] + "</PROJECTION_PARAMETER_VALUE>" + LS,
        /*  7 */ "                        </Projection_Parameter>" + LS,
        /*  8 */ "                        <Projection_Parameter>" + LS,
        /*  9 */ "                            <PROJECTION_PARAMETER_NAME>" + _paramNames[1] + "</PROJECTION_PARAMETER_NAME>" + LS, // notMandatory
        /* 30 */ "                            <PROJECTION_PARAMETER_VALUE unit=\"" + _paramUnit[1] + "\">" + _expValues[1] + "</PROJECTION_PARAMETER_VALUE>" + LS,
        /*  1 */ "                        </Projection_Parameter>" + LS,
        /*  2 */ "                        <Projection_Parameter>" + LS,
        /*  3 */ "                            <PROJECTION_PARAMETER_NAME>" + _paramNames[2] + "</PROJECTION_PARAMETER_NAME>" + LS, // notMandatory
        /*  4 */ "                            <PROJECTION_PARAMETER_VALUE unit=\"" + _paramUnit[2] + "\">" + _expValues[2] + "</PROJECTION_PARAMETER_VALUE>" + LS,
        /*  5 */ "                        </Projection_Parameter>" + LS,
        /*  6 */ "                        <Projection_Parameter>" + LS,
        /*  7 */ "                            <PROJECTION_PARAMETER_NAME>" + _paramNames[3] + "</PROJECTION_PARAMETER_NAME>" + LS, // notMandatory
        /*  8 */ "                            <PROJECTION_PARAMETER_VALUE unit=\"" + _paramUnit[3] + "\">" + _expValues[3] + "</PROJECTION_PARAMETER_VALUE>" + LS,
        /*  9 */ "                        </Projection_Parameter>" + LS,
        /* 40 */ "                        <Projection_Parameter>" + LS,
        /*  1 */ "                            <PROJECTION_PARAMETER_NAME>" + _paramNames[4] + "</PROJECTION_PARAMETER_NAME>" + LS, // notMandatory
        /*  2 */ "                            <PROJECTION_PARAMETER_VALUE unit=\"" + _paramUnit[4] + "\">" + _expValues[4] + "</PROJECTION_PARAMETER_VALUE>" + LS,
        /*  3 */ "                        </Projection_Parameter>" + LS,
        /*  4 */ "                        <Projection_Parameter>" + LS,
        /*  5 */ "                            <PROJECTION_PARAMETER_NAME>" + _paramNames[5] + "</PROJECTION_PARAMETER_NAME>" + LS, // notMandatory
        /*  6 */ "                            <PROJECTION_PARAMETER_VALUE unit=\"" + _paramUnit[5] + "\">" + _expValues[5] + "</PROJECTION_PARAMETER_VALUE>" + LS,
        /*  7 */ "                        </Projection_Parameter>" + LS,
        /*  8 */ "                        <Projection_Parameter>" + LS,
        /*  9 */ "                            <PROJECTION_PARAMETER_NAME>" + _paramNames[6] + "</PROJECTION_PARAMETER_NAME>" + LS, // notMandatory
        /* 50 */ "                            <PROJECTION_PARAMETER_VALUE unit=\"" + _paramUnit[6] + "\">" + _expValues[6] + "</PROJECTION_PARAMETER_VALUE>" + LS,
        /*  1 */ "                        </Projection_Parameter>" + LS,
        /*  2 */ "                    </Projection_Parameters>" + LS,
        /*  3 */ "                </Projection_CT_Method>" + LS,
        /*  4 */ "            </Projection>" + LS,
        /*  5 */ "            <MAP_INFO>" + LS +
        /*  5 */ "                <PIXEL_X value=\"" + _pixelX + "\" />" + LS,
        /*  6 */ "                <PIXEL_Y value=\"" + _pixelY + "\" />" + LS,
        /*  7 */ "                <EASTING value=\"" + _easting + "\" />" + LS,
        /*  8 */ "                <NORTHING value=\"" + _northing + "\" />" + LS,
        /*  9 */ "                <ORIENTATION value=\"" + _orientation + "\" />" + LS,
        /* 60 */ "                <PIXELSIZE_X value=\"" + _pixelSizeX + "\" />" + LS,
        /*  1 */ "                <PIXELSIZE_Y value=\"" + _pixelSizeY + "\" />" + LS,
        /*  2 */ "                <NODATA_VALUE value=\"" + _noDataValue + "\" />" + LS,
        /*  3 */ "                <MAPUNIT value=\"" + _mapUnit + "\" />" + LS,
        /*  4 */ "                <ORTHORECTIFIED value=\"" + _orthorectified + "\" />" + LS,
        /*  5 */ "                <ELEVATION_MODEL value=\"" + _elevModelName + "\" />" + LS,
        /*  6 */ "                <SCENE_FITTED value=\"" + _sceneFitted + "\" />" + LS,
        /*  7 */ "                <SCENE_WIDTH value=\"" + _sceneWidth + "\" />" + LS,
        /*  8 */ "                <SCENE_HEIGHT value=\"" + _sceneHeight + "\" />" + LS,
        /*  9 */ "                <RESAMPLING value=\"" + _resamplingName + "\" />" + LS,
        /* 70 */ "            </MAP_INFO>" + LS,
        /*  1 */ "        </Horizontal_CS>" + LS,
        /*  2 */ "    </Coordinate_Reference_System>" + LS,
        /*  3 */ "</" + DimapProductConstants.TAG_ROOT + ">" + LS
    };

    // fields special for FXYGeoCoding
    private final Datum _datum = Datum.WGS_84;
    private final double[] _xCoefficients = new double[]{0, 1, 2};
    private final double[] _yCoefficients = new double[]{3, 4, 5};
    private final double[] _latCoefficients = new double[]{6, 7, 8};
    private final double[] _lonCoefficients = new double[]{9, 10, 11};
    private final FXYGeoCoding fxyGeoCoding = new FXYGeoCoding(_pixelX, _pixelY, _pixelSizeX, _pixelSizeY,
                                                               FXYSum.createFXYSum(1, _xCoefficients),
                                                               FXYSum.createFXYSum(1, _yCoefficients),
                                                               FXYSum.createFXYSum(1, _latCoefficients),
                                                               FXYSum.createFXYSum(1, _lonCoefficients),
                                                               _datum);

    private final String _xmlFXYGeoCodingString =
            "<" + DimapProductConstants.TAG_ROOT + ">" + LS +
            "    <Coordinate_Reference_System>" + LS +
            "        <Horizontal_CS>" + LS +
            "            <HORIZONTAL_CS_TYPE>GEOGRAPHIC</HORIZONTAL_CS_TYPE>" + LS +
            "            <Geographic_CS>" + LS +
            "                <Horizontal_Datum>" + LS +
            "                    <HORIZONTAL_DATUM_NAME>" + _datum.getName() + "</HORIZONTAL_DATUM_NAME>" + LS +
            "                    <Ellipsoid>" + LS +
            "                        <ELLIPSOID_NAME>" + _datum.getEllipsoid().getName() + "</ELLIPSOID_NAME>" + LS +
            "                        <Ellipsoid_Parameters>" + LS +
            "                            <ELLIPSOID_MAJ_AXIS unit=\"M\">" + _datum.getEllipsoid().getSemiMajor() + "</ELLIPSOID_MAJ_AXIS>" + LS +
            "                            <ELLIPSOID_MIN_AXIS unit=\"M\">" + _datum.getEllipsoid().getSemiMinor() + "</ELLIPSOID_MIN_AXIS>" + LS +
            "                        </Ellipsoid_Parameters>" + LS +
            "                    </Ellipsoid>" + LS +
            "                </Horizontal_Datum>" + LS +
            "            </Geographic_CS>" + LS +
            "        </Horizontal_CS>" + LS +
            "    </Coordinate_Reference_System>" + LS +
            "    <Geoposition>" + LS +
            "        <Geoposition_Insert>" + LS +
            "            <ULXMAP unit=\"M\">" + _pixelX + "</ULXMAP>" + LS +
            "            <ULYMAP unit=\"M\">" + _pixelY + "</ULYMAP>" + LS +
            "            <XDIM unit=\"M\">" + _pixelSizeX + "</XDIM>" + LS +
            "            <YDIM unit=\"M\">" + _pixelSizeY + "</YDIM>" + LS +
            "        </Geoposition_Insert>" + LS +
            "        <Simplified_Location_Model>" + LS +
            "            <Direct_Location_Model order=\"" + fxyGeoCoding.getPixelXFunction().getOrder() + "\">" + LS +
            "                <lc_List>" + LS +
            "                    <lc index=\"0\">" + _lonCoefficients[0] + "</lc>" + LS +
            "                    <lc index=\"1\">" + _lonCoefficients[1] + "</lc>" + LS +
            "                    <lc index=\"2\">" + _lonCoefficients[2] + "</lc>" + LS +
            "                </lc_List>" + LS +
            "                <pc_List>" + LS +
            "                    <pc index=\"0\">" + _latCoefficients[0] + "</pc>" + LS +
            "                    <pc index=\"1\">" + _latCoefficients[1] + "</pc>" + LS +
            "                    <pc index=\"2\">" + _latCoefficients[2] + "</pc>" + LS +
            "                </pc_List>" + LS +
            "            </Direct_Location_Model>" + LS +
            "            <Reverse_Location_Model order=\"" + fxyGeoCoding.getLatFunction().getOrder() + "\">" + LS +
            "                <ic_List>" + LS +
            "                    <ic index=\"0\">" + _xCoefficients[0] + "</ic>" + LS +
            "                    <ic index=\"1\">" + _xCoefficients[1] + "</ic>" + LS +
            "                    <ic index=\"2\">" + _xCoefficients[2] + "</ic>" + LS +
            "                </ic_List>" + LS +
            "                <jc_List>" + LS +
            "                    <jc index=\"0\">" + _yCoefficients[0] + "</jc>" + LS +
            "                    <jc index=\"1\">" + _yCoefficients[1] + "</jc>" + LS +
            "                    <jc index=\"2\">" + _yCoefficients[2] + "</jc>" + LS +
            "                </jc_List>" + LS +
            "            </Reverse_Location_Model>" + LS +
            "        </Simplified_Location_Model>" + LS +
            "    </Geoposition>" + LS +
            "</" + DimapProductConstants.TAG_ROOT + ">";

    private final double[] _xCoefficients1 = new double[]{0, 1, 2};
    private final double[] _yCoefficients1 = new double[]{3, 4, 5};
    private final double[] _latCoefficients1 = new double[]{6, 7, 8};
    private final double[] _lonCoefficients1 = new double[]{9, 10, 11};
    private final FXYGeoCoding _fxyGeoCoding1 = new FXYGeoCoding(_pixelX, _pixelY, _pixelSizeX, _pixelSizeY,
                                                                 FXYSum.createFXYSum(1, _xCoefficients1),
                                                                 FXYSum.createFXYSum(1, _yCoefficients1),
                                                                 FXYSum.createFXYSum(1, _latCoefficients1),
                                                                 FXYSum.createFXYSum(1, _lonCoefficients1),
                                                                 _datum);
    private final double[] _xCoefficients2 = new double[]{12, 13, 14};
    private final double[] _yCoefficients2 = new double[]{15, 16, 17};
    private final double[] _latCoefficients2 = new double[]{18, 19, 20};
    private final double[] _lonCoefficients2 = new double[]{21, 22, 23};
    private final FXYGeoCoding _fxyGeoCoding2 = new FXYGeoCoding(_pixelX, _pixelY, _pixelSizeX, _pixelSizeY,
                                                                 FXYSum.createFXYSum(1, _xCoefficients2),
                                                                 FXYSum.createFXYSum(1, _yCoefficients2),
                                                                 FXYSum.createFXYSum(1, _latCoefficients2),
                                                                 FXYSum.createFXYSum(1, _lonCoefficients2),
                                                                 _datum);

    private final String _xmlBandedFXYGeoCodingString =
            "<" + DimapProductConstants.TAG_ROOT + ">" + LS +
            "    <Coordinate_Reference_System>" + LS +
            "        <Horizontal_CS>" + LS +
            "            <HORIZONTAL_CS_TYPE>GEOGRAPHIC</HORIZONTAL_CS_TYPE>" + LS +
            "            <Geographic_CS>" + LS +
            "                <Horizontal_Datum>" + LS +
            "                    <HORIZONTAL_DATUM_NAME>" + _datum.getName() + "</HORIZONTAL_DATUM_NAME>" + LS +
            "                    <Ellipsoid>" + LS +
            "                        <ELLIPSOID_NAME>" + _datum.getEllipsoid().getName() + "</ELLIPSOID_NAME>" + LS +
            "                        <Ellipsoid_Parameters>" + LS +
            "                            <ELLIPSOID_MAJ_AXIS unit=\"M\">" + _datum.getEllipsoid().getSemiMajor() + "</ELLIPSOID_MAJ_AXIS>" + LS +
            "                            <ELLIPSOID_MIN_AXIS unit=\"M\">" + _datum.getEllipsoid().getSemiMinor() + "</ELLIPSOID_MIN_AXIS>" + LS +
            "                        </Ellipsoid_Parameters>" + LS +
            "                    </Ellipsoid>" + LS +
            "                </Horizontal_Datum>" + LS +
            "            </Geographic_CS>" + LS +
            "        </Horizontal_CS>" + LS +
            "    </Coordinate_Reference_System>" + LS +
            "    <Geoposition>" + LS +
            "        <BAND_INDEX>0</BAND_INDEX>" + LS +
            "        <Geoposition_Insert>" + LS +
            "             <ULXMAP>" + _pixelX + "</ULXMAP>" + LS +
            "             <ULYMAP>" + _pixelY + "</ULYMAP>" + LS +
            "             <XDIM>" + _pixelSizeX + "</XDIM>" + LS +
            "             <YDIM>" + _pixelSizeY + "</YDIM>" + LS +
            "        </Geoposition_Insert>" + LS +
            "         <Simplified_Location_Model>" + LS +
            "             <Direct_Location_Model order=\"" + _fxyGeoCoding1.getPixelXFunction().getOrder() + "\">" + LS +
            "                <lc_List>" + LS +
            "                    <lc index=\"0\">" + _lonCoefficients1[0] + "</lc>" + LS +
            "                    <lc index=\"1\">" + _lonCoefficients1[1] + "</lc>" + LS +
            "                    <lc index=\"2\">" + _lonCoefficients1[2] + "</lc>" + LS +
            "                </lc_List>" + LS +
            "                <pc_List>" + LS +
            "                    <pc index=\"0\">" + _latCoefficients1[0] + "</pc>" + LS +
            "                    <pc index=\"1\">" + _latCoefficients1[1] + "</pc>" + LS +
            "                    <pc index=\"2\">" + _latCoefficients1[2] + "</pc>" + LS +
            "                </pc_List>" + LS +
            "            </Direct_Location_Model>" + LS +
            "            <Reverse_Location_Model order=\"" + _fxyGeoCoding1.getLatFunction().getOrder() + "\">" + LS +
            "                <ic_List>" + LS +
            "                    <ic index=\"0\">" + _xCoefficients1[0] + "</ic>" + LS +
            "                    <ic index=\"1\">" + _xCoefficients1[1] + "</ic>" + LS +
            "                    <ic index=\"2\">" + _xCoefficients1[2] + "</ic>" + LS +
            "                </ic_List>" + LS +
            "                <jc_List>" + LS +
            "                    <jc index=\"0\">" + _yCoefficients1[0] + "</jc>" + LS +
            "                    <jc index=\"1\">" + _yCoefficients1[1] + "</jc>" + LS +
            "                    <jc index=\"2\">" + _yCoefficients1[2] + "</jc>" + LS +
            "                </jc_List>" + LS +
            "            </Reverse_Location_Model>" + LS +
            "        </Simplified_Location_Model>" + LS +
            "    </Geoposition>" + LS +
            "    <Geoposition>" + LS +
            "        <BAND_INDEX>1</BAND_INDEX>" + LS +
            "        <Geoposition_Insert>" + LS +
            "             <ULXMAP unit=\"M\">" + _pixelX + "</ULXMAP>" + LS +
            "             <ULYMAP unit=\"M\">" + _pixelY + "</ULYMAP>" + LS +
            "             <XDIM unit=\"M\">" + _pixelSizeX + "</XDIM>" + LS +
            "             <YDIM unit=\"M\">" + _pixelSizeY + "</YDIM>" + LS +
            "        </Geoposition_Insert>" + LS +
            "         <Simplified_Location_Model>" + LS +
            "             <Direct_Location_Model order=\"" + _fxyGeoCoding2.getPixelXFunction().getOrder() + "\">" + LS +
            "                <lc_List>" + LS +
            "                    <lc index=\"0\">" + _lonCoefficients2[0] + "</lc>" + LS +
            "                    <lc index=\"1\">" + _lonCoefficients2[1] + "</lc>" + LS +
            "                    <lc index=\"2\">" + _lonCoefficients2[2] + "</lc>" + LS +
            "                </lc_List>" + LS +
            "                <pc_List>" + LS +
            "                    <pc index=\"0\">" + _latCoefficients2[0] + "</pc>" + LS +
            "                    <pc index=\"1\">" + _latCoefficients2[1] + "</pc>" + LS +
            "                    <pc index=\"2\">" + _latCoefficients2[2] + "</pc>" + LS +
            "                </pc_List>" + LS +
            "            </Direct_Location_Model>" + LS +
            "            <Reverse_Location_Model order=\"" + _fxyGeoCoding2.getLatFunction().getOrder() + "\">" + LS +
            "                <ic_List>" + LS +
            "                    <ic index=\"0\">" + _xCoefficients2[0] + "</ic>" + LS +
            "                    <ic index=\"1\">" + _xCoefficients2[1] + "</ic>" + LS +
            "                    <ic index=\"2\">" + _xCoefficients2[2] + "</ic>" + LS +
            "                </ic_List>" + LS +
            "                <jc_List>" + LS +
            "                    <jc index=\"0\">" + _yCoefficients2[0] + "</jc>" + LS +
            "                    <jc index=\"1\">" + _yCoefficients2[1] + "</jc>" + LS +
            "                    <jc index=\"2\">" + _yCoefficients2[2] + "</jc>" + LS +
            "                </jc_List>" + LS +
            "            </Reverse_Location_Model>" + LS +
            "        </Simplified_Location_Model>" + LS +
            "    </Geoposition>" + LS +
            "</" + DimapProductConstants.TAG_ROOT + ">";


    private Product product;

    @Override
    protected void setUp() throws Exception {
        product = new Product("product", "type", 200, 300);
        product.addBand("b1", ProductData.TYPE_INT8);
        product.addBand("b2", ProductData.TYPE_INT8);
    }

    public void testCreateGeoCodingForMapProjectionWithOldDimapFormat() {
        final String projectionName = "Lambert Conformal Conic";

        final double pixelX = 1.0;
        final double pixelY = 2.0;
        final double easting = 3.0;
        final double northing = 4.0;
        final double pixelSizeX = 5.0;
        final double pixelSizeY = 6.0;
        final String datumName = "WGS-84";
        final int sceneWidth = 140;
        final int sceneHeight = 232;

        final String xmlMapGeocodingStringStyleOldFormat =
                "<" + DimapProductConstants.TAG_ROOT + ">" + LS +
                /*  1 */ "    <Coordinate_Reference_System>" + LS +
                /*  2 */ "        <Geocoding_Map>" + LS +
                /*  3 */ "            <MAP_INFO>" + projectionName + ", " + pixelX + ", " + pixelY + ", " + easting + ", " + northing + ", " + pixelSizeX + ", " + pixelSizeY + ", " + datumName + ", units=meter, " + sceneWidth + ", " + sceneHeight + "</MAP_INFO>" + LS +
                /*  4 */ "        </Geocoding_Map>" + LS +
                /*  5 */ "    </Coordinate_Reference_System>" + LS +
                /*  6 */ "</" + DimapProductConstants.TAG_ROOT + ">" + LS;

        final Ellipsoid expEllipsoid = Ellipsoid.WGS_84;
        final LambertConformalConicDescriptor expDescriptor = new LambertConformalConicDescriptor();
        final MapTransform expTransform = expDescriptor.createTransform(null);
        final Document dom = DimapProductHelpers.createDom(
                new ByteArrayInputStream(xmlMapGeocodingStringStyleOldFormat.getBytes()));

        //test
        GeoCoding[] geoCodings = DimapProductHelpers.createGeoCoding(dom, product);
        assertNotNull(geoCodings);
        assertEquals(1, geoCodings.length);

        final GeoCoding geoCoding = geoCodings[0];

        assertNotNull(geoCoding);
        assertNotNull(geoCoding);
        assertEquals(MapGeoCoding.class, geoCoding.getClass());
        final MapGeoCoding mapGeoCoding = (MapGeoCoding) geoCoding;

        final MapInfo actualMapinfo = mapGeoCoding.getMapInfo();
        assertEquals(pixelX, actualMapinfo.getPixelX(), 1.0e-6);
        assertEquals(pixelY, actualMapinfo.getPixelY(), 1.0e-6);
        assertEquals(easting, actualMapinfo.getEasting(), 1.0e-6);
        assertEquals(northing, actualMapinfo.getNorthing(), 1.0e-6);
        assertEquals(pixelSizeX, actualMapinfo.getPixelSizeX(), 1.0e-6);
        assertEquals(pixelSizeY, actualMapinfo.getPixelSizeY(), 1.0e-6);
        assertEquals(sceneWidth, actualMapinfo.getSceneWidth());
        assertEquals(sceneHeight, actualMapinfo.getSceneHeight());

        final Datum actualDatum = actualMapinfo.getDatum();
        assertEquals(datumName, actualDatum.getName());

        final MapProjection actualProjection = actualMapinfo.getMapProjection();
        assertEquals(projectionName, actualProjection.getName());

        final Ellipsoid actualEllipsoid = actualDatum.getEllipsoid();
        assertEquals(expEllipsoid.getName(), actualEllipsoid.getName());
        assertEquals(expEllipsoid.getSemiMajor(), actualEllipsoid.getSemiMajor(), 1.0e-10);
        assertEquals(expEllipsoid.getSemiMinor(), actualEllipsoid.getSemiMinor(), 1.0e-10);

        final MapTransform actualTransform = actualProjection.getMapTransform();
        assertEquals(expTransform.getClass(), actualTransform.getClass());

        final double[] actualValues = actualTransform.getParameterValues();
        assertEquals(true, Arrays.equals(expTransform.getParameterValues(), actualValues));

        final MapTransformDescriptor actualDescriptor = actualTransform.getDescriptor();
        assertEquals(expDescriptor.getClass(), actualDescriptor.getClass());
    }

    public void testCreateGeoCodingForMapProjectionWithFullValidDimap_1_4_0_Format() {
        final int[] allLines = new int[0];

        final String fullXMLString = createXMLString(allLines);
        final Document dom = DimapProductHelpers.createDom(new ByteArrayInputStream(fullXMLString.getBytes()));
        final GeoCoding geoCoding = DimapProductHelpers.createGeoCoding(dom, product)[0];

        assertNotNull(geoCoding);
        assertEquals(MapGeoCoding.class, geoCoding.getClass());
        final MapGeoCoding mapGeoCoding = (MapGeoCoding) geoCoding;

        final MapInfo actualMapinfo = mapGeoCoding.getMapInfo();
        assertEquals(_pixelX, actualMapinfo.getPixelX(), 1.0e-6);
        assertEquals(_pixelY, actualMapinfo.getPixelY(), 1.0e-6);
        assertEquals(_easting, actualMapinfo.getEasting(), 1.0e-6);
        assertEquals(_northing, actualMapinfo.getNorthing(), 1.0e-6);
        assertEquals(_pixelSizeX, actualMapinfo.getPixelSizeX(), 1.0e-6);
        assertEquals(_pixelSizeY, actualMapinfo.getPixelSizeY(), 1.0e-6);
        assertEquals(_sceneWidth, actualMapinfo.getSceneWidth());
        assertEquals(_sceneHeight, actualMapinfo.getSceneHeight());

        final Datum actualDatum = actualMapinfo.getDatum();
        assertEquals(_datumName, actualDatum.getName());

        final MapProjection actualProjection = actualMapinfo.getMapProjection();
        assertEquals(_projectionName, actualProjection.getName());
        assertEquals(_mapUnit, actualProjection.getMapUnit());

        final Ellipsoid actualEllipsoid = actualDatum.getEllipsoid();
        assertEquals(_ellipsoidName, actualEllipsoid.getName());
        assertEquals(_semiMajor, actualEllipsoid.getSemiMajor(), 1.0e-10);
        assertEquals(_semiMinor, actualEllipsoid.getSemiMinor(), 1.0e-10);

        final MapTransform actualTransform = actualProjection.getMapTransform();
        assertEquals(LambertConformalConicDescriptor.LCCT.class, actualTransform.getClass());

        final double[] actualValues = actualTransform.getParameterValues();
        for (int i = 0; i < actualValues.length; i++) {
            assertEquals(_expValues[i], actualValues[i], 1.0e-5);
        }

        final MapTransformDescriptor actualDescriptor = actualTransform.getDescriptor();
        assertEquals(LambertConformalConicDescriptor.class, actualDescriptor.getClass());
    }

    public void testCreateGeoCodingForMapProjectionWithValidDimap_1_4_0_Format() {

        final String xmlStringWithoutNotMandatoryLines = createXMLString(_notMandatoryLines);
        final byte[] bytes = xmlStringWithoutNotMandatoryLines.getBytes();
        final Document dom = DimapProductHelpers.createDom(new ByteArrayInputStream(bytes));
        final GeoCoding geoCoding = DimapProductHelpers.createGeoCoding(dom, product)[0];

        assertNotNull(geoCoding);
        assertEquals(MapGeoCoding.class, geoCoding.getClass());
        final MapGeoCoding mapGeoCoding = (MapGeoCoding) geoCoding;

        final MapInfo actualMapinfo = mapGeoCoding.getMapInfo();
        assertEquals(_pixelX, actualMapinfo.getPixelX(), 1.0e-6);
        assertEquals(_pixelY, actualMapinfo.getPixelY(), 1.0e-6);
        assertEquals(_easting, actualMapinfo.getEasting(), 1.0e-6);
        assertEquals(_northing, actualMapinfo.getNorthing(), 1.0e-6);
        assertEquals(_pixelSizeX, actualMapinfo.getPixelSizeX(), 1.0e-6);
        assertEquals(_pixelSizeY, actualMapinfo.getPixelSizeY(), 1.0e-6);
        assertEquals(_sceneWidth, actualMapinfo.getSceneWidth());
        assertEquals(_sceneHeight, actualMapinfo.getSceneHeight());

        final Datum actualDatum = actualMapinfo.getDatum();
        assertEquals(_datumName, actualDatum.getName());

        final MapProjection actualProjection = actualMapinfo.getMapProjection();
        assertEquals(_projectionName, actualProjection.getName());
        assertEquals(_mapUnit, actualProjection.getMapUnit());

        final Ellipsoid actualEllipsoid = actualDatum.getEllipsoid();
        assertEquals(_ellipsoidName, actualEllipsoid.getName());
        assertEquals(_semiMajor, actualEllipsoid.getSemiMajor(), 1.0e-10);
        assertEquals(_semiMinor, actualEllipsoid.getSemiMinor(), 1.0e-10);

        final MapTransform actualTransform = actualProjection.getMapTransform();
        assertEquals(LambertConformalConicDescriptor.LCCT.class, actualTransform.getClass());

        final double[] actualValues = actualTransform.getParameterValues();
        assertEquals(true, Arrays.equals(_expValues, actualValues));

        final MapTransformDescriptor actualDescriptor = actualTransform.getDescriptor();
        assertEquals(LambertConformalConicDescriptor.class, actualDescriptor.getClass());
    }

    public void testCreateGeoCodingForMapProjectionWithInvalidDimap_1_4_0_Format() {
        Document dom;

        // without ellipsoid name Element
        dom = createDom(withoutMandatoryLines(new int[]{11}));
        assertEquals(null, DimapProductHelpers.createGeoCoding(dom, product));

        // without ellipsoid major axis Element
        dom = createDom(withoutMandatoryLines(new int[]{13}));
        assertEquals(null, DimapProductHelpers.createGeoCoding(dom, product));

        // without ellipsoid minor axis Element
        dom = createDom(withoutMandatoryLines(new int[]{14}));
        assertEquals(null, DimapProductHelpers.createGeoCoding(dom, product));

        // without ellipsoid parameters Element
        dom = createDom(withoutMandatoryLines(ArrayUtils.createIntArray(12, 15)));
        assertEquals(null, DimapProductHelpers.createGeoCoding(dom, product));

        // without ellipsoid Element
        dom = createDom(withoutMandatoryLines(ArrayUtils.createIntArray(10, 16)));
        assertEquals(null, DimapProductHelpers.createGeoCoding(dom, product));

        // without horizontal datum name Element
        dom = createDom(withoutMandatoryLines(new int[]{9}));
        assertEquals(null, DimapProductHelpers.createGeoCoding(dom, product));

        // without horizontal datum Element
        dom = createDom(withoutMandatoryLines(ArrayUtils.createIntArray(8, 17)));
        assertEquals(null, DimapProductHelpers.createGeoCoding(dom, product));

        // without Geographic_CS Element
        dom = createDom(withoutMandatoryLines(ArrayUtils.createIntArray(6, 18)));
        assertEquals(null, DimapProductHelpers.createGeoCoding(dom, product));

        // without Horizontal_CS Element
        dom = createDom(withoutMandatoryLines(ArrayUtils.createIntArray(3, 71)));
        assertEquals(null, DimapProductHelpers.createGeoCoding(dom, product));

        // without Coordinate_Reference_System Element
        dom = createDom(withoutMandatoryLines(ArrayUtils.createIntArray(1, 72)));
        assertEquals(null, DimapProductHelpers.createGeoCoding(dom, product));

        // without Projection Element
        dom = createDom(withoutMandatoryLines(ArrayUtils.createIntArray(19, 54)));
        assertEquals(null, DimapProductHelpers.createGeoCoding(dom, product));

        // without NAME Element
        dom = createDom(withoutMandatoryLines(new int[]{20}));
        assertEquals(null, DimapProductHelpers.createGeoCoding(dom, product));

        // without Projection_CT_Method Element
        dom = createDom(withoutMandatoryLines(ArrayUtils.createIntArray(21, 53)));
        assertEquals(null, DimapProductHelpers.createGeoCoding(dom, product));

        // without PROJECTION_CT_NAME Element
        dom = createDom(withoutMandatoryLines(new int[]{22}));
        assertEquals(null, DimapProductHelpers.createGeoCoding(dom, product));

        // without Projection_Parameters Element
        dom = createDom(withoutMandatoryLines(ArrayUtils.createIntArray(23, 52)));
        assertEquals(null, DimapProductHelpers.createGeoCoding(dom, product));

        // without all Projection_Parameter Elements
        dom = createDom(withoutMandatoryLines(ArrayUtils.createIntArray(24, 51)));
        assertEquals(null, DimapProductHelpers.createGeoCoding(dom, product));

        // without PROJECTION_PARAMETER_VALUE Element
        dom = createDom(withoutMandatoryLines(new int[]{26}));
        assertEquals(null, DimapProductHelpers.createGeoCoding(dom, product));

        // without MAP_INFO Element
        dom = createDom(withoutMandatoryLines(ArrayUtils.createIntArray(55, 70)));
        assertEquals(null, DimapProductHelpers.createGeoCoding(dom, product));
    }

    public void testCreateGeoCodingForFXYGeoCoding() {
        final byte[] bytes = _xmlFXYGeoCodingString.getBytes();
        final Document dom = DimapProductHelpers.createDom(new ByteArrayInputStream(bytes));
        final GeoCoding geoCoding = DimapProductHelpers.createGeoCoding(dom, product)[0];

        assertNotNull(geoCoding);
        assertEquals(FXYGeoCoding.class, geoCoding.getClass());
        assertEqual(this.fxyGeoCoding, (FXYGeoCoding) geoCoding);

    }

    public void testCreateGeoCodingForPixelGeoCoding() throws IOException {
        final Band latBand = product.getBand("b1");
        final Band lonBand = product.getBand("b2");
        final byte[] bandData = new byte[product.getSceneRasterWidth() * product.getSceneRasterHeight()];
        latBand.setDataElems(bandData);
        lonBand.setDataElems(bandData);
        final PixelGeoCoding sourcePixelGeoCoding = new PixelGeoCoding(latBand, lonBand,
                                                                       "Not NaN", 4, ProgressMonitor.NULL);
        final byte[] bytes = createPixelGeoCodingString(sourcePixelGeoCoding).getBytes();
        final Document dom = DimapProductHelpers.createDom(new ByteArrayInputStream(bytes));
        final GeoCoding geoCoding = DimapProductHelpers.createGeoCoding(dom, product)[0];

        assertNotNull(geoCoding);
        assertEquals(PixelGeoCoding.class, geoCoding.getClass());
        final PixelGeoCoding pixelGeoCoding = (PixelGeoCoding) geoCoding;

        assertEquals("b1", pixelGeoCoding.getLatBand().getName());
        assertEquals("b2", pixelGeoCoding.getLonBand().getName());
        assertEquals("Not NaN", pixelGeoCoding.getValidMask());
        assertEquals(4, pixelGeoCoding.getSearchRadius());
        assertNotNull(pixelGeoCoding.getPixelPosEstimator());
        assertEqual(fxyGeoCoding, (FXYGeoCoding) pixelGeoCoding.getPixelPosEstimator());

    }

    public void testCreateGeoCodingForCrsGeoCoding() throws Exception {
        final Rectangle2D.Double imageBounds = new Rectangle2D.Double(0, 0, product.getSceneRasterWidth(),
                                                                      product.getSceneRasterHeight());
        final AffineTransform expectedI2m = new AffineTransform(0.12, 1.23, 2.34, 3.45, 4.56, 5.67);
        final CoordinateReferenceSystem expectedCrs = CRS.decode("EPSG:4326");
        final byte[] bytes = createCrsGeoCodingString(new CrsGeoCoding(expectedCrs, imageBounds, expectedI2m)).getBytes();
        final Document dom = DimapProductHelpers.createDom(new ByteArrayInputStream(bytes));
        final GeoCoding geoCoding = DimapProductHelpers.createGeoCoding(dom, product)[0];

        assertNotNull(geoCoding);
        assertEquals(CrsGeoCoding.class, geoCoding.getClass());
        final CrsGeoCoding crsGeoCoding = (CrsGeoCoding) geoCoding;

        final CoordinateReferenceSystem modelCRS = crsGeoCoding.getModelCRS();
        // ignoring metadata because scope and domainOfValidity are not restored
        // but not important for our GeoCoding
        assertTrue(CRS.equalsIgnoreMetadata(expectedCrs, modelCRS));
        assertEquals(expectedI2m, crsGeoCoding.getImageToModelTransform());
    }

    public void testReadingGeoCodingPerBand() {
        final byte[] bytes = _xmlBandedFXYGeoCodingString.getBytes();
        final Document dom = DimapProductHelpers.createDom(new ByteArrayInputStream(bytes));
        final GeoCoding[] geoCodings = DimapProductHelpers.createGeoCoding(dom, product);

        assertFalse(product.isUsingSingleGeoCoding());
        assertTrue(product.getGeoCoding() == null);

        final GeoCoding geoCoding1 = geoCodings[0];
        final GeoCoding geoCoding2 = geoCodings[1];
        assertTrue(geoCoding1 instanceof FXYGeoCoding);
        assertTrue(geoCoding2 instanceof FXYGeoCoding);
        assertNotSame(geoCoding1, geoCoding2);

        final FXYGeoCoding fxyGeoCoding1 = (FXYGeoCoding) geoCoding1;
        final FXYGeoCoding fxyGeoCoding2 = (FXYGeoCoding) geoCoding2;

        assertEqual(_fxyGeoCoding1, fxyGeoCoding1);
        assertEqual(_fxyGeoCoding2, fxyGeoCoding2);
    }

    private void assertEqual(final FXYGeoCoding expectedGeoCoding, final FXYGeoCoding actualGeoCoding) {
        assertEquals(expectedGeoCoding.getDatum().getName(), actualGeoCoding.getDatum().getName());
        assertEquals(expectedGeoCoding.getDatum().getEllipsoid().getName(),
                     actualGeoCoding.getDatum().getEllipsoid().getName());
        assertEquals(expectedGeoCoding.getDatum().getEllipsoid().getSemiMajor(),
                     actualGeoCoding.getDatum().getEllipsoid().getSemiMajor(), 1.0e-6);
        assertEquals(expectedGeoCoding.getDatum().getEllipsoid().getSemiMinor(),
                     actualGeoCoding.getDatum().getEllipsoid().getSemiMinor(), 1.0e-6);
        assertEquals(expectedGeoCoding.getLatFunction().getOrder(), actualGeoCoding.getLatFunction().getOrder());
        assertTrue(ArrayUtils.equalArrays(expectedGeoCoding.getLatFunction().getCoefficients(),
                                          actualGeoCoding.getLatFunction().getCoefficients(), 1.0e-6));
        assertEquals(expectedGeoCoding.getLonFunction().getOrder(), actualGeoCoding.getLonFunction().getOrder());
        assertTrue(ArrayUtils.equalArrays(expectedGeoCoding.getLonFunction().getCoefficients(),
                                          actualGeoCoding.getLonFunction().getCoefficients(), 1.0e-6));
        assertEquals(expectedGeoCoding.getPixelXFunction().getOrder(), actualGeoCoding.getPixelXFunction().getOrder());
        assertTrue(ArrayUtils.equalArrays(expectedGeoCoding.getPixelXFunction().getCoefficients(),
                                          actualGeoCoding.getPixelXFunction().getCoefficients(), 1.0e-6));
        assertEquals(expectedGeoCoding.getPixelYFunction().getOrder(), actualGeoCoding.getPixelYFunction().getOrder());
        assertTrue(ArrayUtils.equalArrays(expectedGeoCoding.getPixelYFunction().getCoefficients(),
                                          actualGeoCoding.getPixelYFunction().getCoefficients(), 1.0e-6));
        assertEquals(expectedGeoCoding.getPixelOffsetX(), actualGeoCoding.getPixelOffsetX(), 1.0e-6);
        assertEquals(expectedGeoCoding.getPixelOffsetY(), actualGeoCoding.getPixelOffsetY(), 1.0e-6);
        assertEquals(expectedGeoCoding.getPixelSizeX(), actualGeoCoding.getPixelSizeX(), 1.0e-6);
        assertEquals(expectedGeoCoding.getPixelSizeY(), actualGeoCoding.getPixelSizeY(), 1.0e-6);
    }

    public void testConvertBeamUnitToDimapUnit() {
        assertEquals("M", DimapProductHelpers.convertBeamUnitToDimapUnit("meter"));
        assertEquals("M", DimapProductHelpers.convertBeamUnitToDimapUnit("meters"));

        assertEquals("KM", DimapProductHelpers.convertBeamUnitToDimapUnit("kilometer"));
        assertEquals("KM", DimapProductHelpers.convertBeamUnitToDimapUnit("kilometers"));

        assertEquals("DEG", DimapProductHelpers.convertBeamUnitToDimapUnit("deg"));
        assertEquals("DEG", DimapProductHelpers.convertBeamUnitToDimapUnit("DEG"));
        assertEquals("DEG", DimapProductHelpers.convertBeamUnitToDimapUnit("degree"));

        assertEquals("RAD", DimapProductHelpers.convertBeamUnitToDimapUnit("rad"));
        assertEquals("RAD", DimapProductHelpers.convertBeamUnitToDimapUnit("RADIAN"));
        assertEquals("RAD", DimapProductHelpers.convertBeamUnitToDimapUnit("radian"));

    }

    public void testConvertDimaUnitToBeamUnit() {
        assertEquals("meter", DimapProductHelpers.convertDimapUnitToBeamUnit("M"));

        assertEquals("kilometer", DimapProductHelpers.convertDimapUnitToBeamUnit("KM"));

        assertEquals("deg", DimapProductHelpers.convertDimapUnitToBeamUnit("DEG"));

        assertEquals("rad", DimapProductHelpers.convertDimapUnitToBeamUnit("RAD"));

    }

    private int[] withoutMandatoryLines(final int[] value) {
        return ArrayUtils.addArrays(_notMandatoryLines, value);
    }

    private Document createDom(final int[] withoutMandatoryLine) {
        final String xmlStringWithoutNotMandatoryLines = createXMLString(withoutMandatoryLine);
        final byte[] bytes = xmlStringWithoutNotMandatoryLines.getBytes();
        return DimapProductHelpers.createDom(new ByteArrayInputStream(bytes));
    }

    private String createXMLString(final int[] notAccepted) {
        Arrays.sort(notAccepted);
        final Vector<String> vector = new Vector<String>();
        for (int i = 0; i < _xmlMapGeocodingStringStyleV1_4_0.length; i++) {
            if (Arrays.binarySearch(notAccepted, i) < 0) {
                vector.add(_xmlMapGeocodingStringStyleV1_4_0[i]);
            }
        }
        final String[] strings = vector.toArray(new String[vector.size()]);
        return StringUtils.arrayToString(strings, "");
    }

    private String createPixelGeoCodingString(PixelGeoCoding geoCoding) {
        return "<" + DimapProductConstants.TAG_ROOT + ">" + LS +
               "    <Geoposition>" + LS +
               "        <LATITUDE_BAND>" + geoCoding.getLatBand().getName() + "</LATITUDE_BAND>" + LS +
               "        <LONGITUDE_BAND>" + geoCoding.getLonBand().getName() + "</LONGITUDE_BAND>" + LS +
               "        <VALID_MASK_EXPRESSION>" + geoCoding.getValidMask() + "</VALID_MASK_EXPRESSION>" + LS +
               "        <SEARCH_RADIUS>" + geoCoding.getSearchRadius() + "</SEARCH_RADIUS>" + LS +
               "        <Pixel_Position_Estimator>" + LS +
               "            <Coordinate_Reference_System>" + LS +
               "                <Horizontal_CS>" + LS +
               "                    <HORIZONTAL_CS_TYPE>GEOGRAPHIC</HORIZONTAL_CS_TYPE>" + LS +
               "                    <Geographic_CS>" + LS +
               "                        <Horizontal_Datum>" + LS +
               "                            <HORIZONTAL_DATUM_NAME>" + _datum.getName() + "</HORIZONTAL_DATUM_NAME>" + LS +
               "                            <Ellipsoid>" + LS +
               "                                <ELLIPSOID_NAME>" + _datum.getEllipsoid().getName() + "</ELLIPSOID_NAME>" + LS +
               "                                <Ellipsoid_Parameters>" + LS +
               "                                    <ELLIPSOID_MAJ_AXIS unit=\"M\">" + _datum.getEllipsoid().getSemiMajor() + "</ELLIPSOID_MAJ_AXIS>" + LS +
               "                                    <ELLIPSOID_MIN_AXIS unit=\"M\">" + _datum.getEllipsoid().getSemiMinor() + "</ELLIPSOID_MIN_AXIS>" + LS +
               "                                </Ellipsoid_Parameters>" + LS +
               "                            </Ellipsoid>" + LS +
               "                        </Horizontal_Datum>" + LS +
               "                    </Geographic_CS>" + LS +
               "                </Horizontal_CS>" + LS +
               "            </Coordinate_Reference_System>" + LS +
               "            <Geoposition>" + LS +
               "                <Geoposition_Insert>" + LS +
               "                    <ULXMAP unit=\"M\">" + _pixelX + "</ULXMAP>" + LS +
               "                    <ULYMAP unit=\"M\">" + _pixelY + "</ULYMAP>" + LS +
               "                    <XDIM unit=\"M\">" + _pixelSizeX + "</XDIM>" + LS +
               "                    <YDIM unit=\"M\">" + _pixelSizeY + "</YDIM>" + LS +
               "                </Geoposition_Insert>" + LS +
               "                <Simplified_Location_Model>" + LS +
               "                    <Direct_Location_Model order=\"" + fxyGeoCoding.getPixelXFunction().getOrder() + "\">" + LS +
               "                        <lc_List>" + LS +
               "                            <lc index=\"0\">" + _lonCoefficients[0] + "</lc>" + LS +
               "                            <lc index=\"1\">" + _lonCoefficients[1] + "</lc>" + LS +
               "                            <lc index=\"2\">" + _lonCoefficients[2] + "</lc>" + LS +
               "                        </lc_List>" + LS +
               "                        <pc_List>" + LS +
               "                            <pc index=\"0\">" + _latCoefficients[0] + "</pc>" + LS +
               "                            <pc index=\"1\">" + _latCoefficients[1] + "</pc>" + LS +
               "                            <pc index=\"2\">" + _latCoefficients[2] + "</pc>" + LS +
               "                        </pc_List>" + LS +
               "                    </Direct_Location_Model>" + LS +
               "                    <Reverse_Location_Model order=\"" + fxyGeoCoding.getLatFunction().getOrder() + "\">" + LS +
               "                        <ic_List>" + LS +
               "                            <ic index=\"0\">" + _xCoefficients[0] + "</ic>" + LS +
               "                            <ic index=\"1\">" + _xCoefficients[1] + "</ic>" + LS +
               "                            <ic index=\"2\">" + _xCoefficients[2] + "</ic>" + LS +
               "                        </ic_List>" + LS +
               "                        <jc_List>" + LS +
               "                            <jc index=\"0\">" + _yCoefficients[0] + "</jc>" + LS +
               "                            <jc index=\"1\">" + _yCoefficients[1] + "</jc>" + LS +
               "                            <jc index=\"2\">" + _yCoefficients[2] + "</jc>" + LS +
               "                        </jc_List>" + LS +
               "                    </Reverse_Location_Model>" + LS +
               "                </Simplified_Location_Model>" + LS +
               "            </Geoposition>" + LS +
               "        </Pixel_Position_Estimator>" + LS +
               "    </Geoposition>" + LS +
               "</" + DimapProductConstants.TAG_ROOT + ">";
    }

    private String createCrsGeoCodingString(CrsGeoCoding geoCoding) {
        final double[] matrix = new double[6];
        geoCoding.getImageToModelTransform().getMatrix(matrix);

        return "<" + DimapProductConstants.TAG_ROOT + ">" + LS +
               "    <Coordinate_Reference_System>" + LS +
               "        <WKT>" + LS +
               geoCoding.getModelCRS().toString() +
               "        </WKT>" + LS +
               "    </Coordinate_Reference_System>" + LS +
               "    <Geoposition>" + LS +
               "        <IMAGE_TO_MODEL_TRANSFORM>"+StringUtils.arrayToCsv(matrix)+"</IMAGE_TO_MODEL_TRANSFORM>" + LS +
               "    </Geoposition>" + LS +
               "</" + DimapProductConstants.TAG_ROOT + ">";
    }

}
