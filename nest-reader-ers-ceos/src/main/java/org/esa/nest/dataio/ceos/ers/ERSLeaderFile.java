package org.esa.nest.dataio.ceos.ers;

import org.esa.nest.dataio.ceos.CeosFileReader;
import org.esa.nest.dataio.ceos.IllegalCeosFormatException;
import org.esa.nest.dataio.ceos.records.BaseRecord;
import org.esa.nest.dataio.ceos.records.BaseSceneHeaderRecord;
import org.esa.beam.framework.datamodel.MetadataAttribute;
import org.esa.beam.framework.datamodel.MetadataElement;
import org.esa.beam.framework.datamodel.ProductData;

import javax.imageio.stream.ImageInputStream;
import java.io.IOException;
import java.util.Calendar;

/*
 * $Id: ERSLeaderFile.java,v 1.5 2008-06-19 17:42:37 lveci Exp $
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

/**
 * This class represents a leader file of an Avnir-2 product.
 *
 * @author Marco Peters
 * @version $Revision: 1.5 $ $Date: 2008-06-19 17:42:37 $
 */
class ERSLeaderFile {

    private static final String UNIT_METER = "meter";
    private static final String UNIT_KILOMETER = "kilometer";
    private static final String UNIT_DEGREE = "degree";
    private static final String UNIT_SECOND = "second";
    private static final String UNIT_METER_PER_SECOND = "m/sec";
    private static final String UNIT_DEGREE_PER_SECOND = "deg/sec";

    private static final String PROJECTION_KEY_RAW = "NNNNN";
    private static final String PROJECTION_KEY_UTM = "YNNNN";
    private static final String PROJECTION_KEY_PS = "NNNNY";

    public final BaseRecord _leaderFDR;
    public final BaseSceneHeaderRecord _sceneHeaderRecord;
    public final BaseRecord _mapProjRecord;
    public final BaseRecord _platformPositionRecord;
    public final BaseRecord _facilityRecord;
    public CeosFileReader _reader;

    private static String mission = "ers";
    private static String leader_recordDefinitionFile = "leader_file.xml";
    private static String scene_recordDefinitionFile = "scene_record.xml";
    private static String mapproj_recordDefinitionFile = "map_proj_record.xml";
    private static String platform_recordDefinitionFile = "platform_position_record.xml";
    private static String facility_recordDefinitionFile = "facility_record.xml";

    public ERSLeaderFile(final ImageInputStream leaderStream) throws IOException,
                                                                        IllegalCeosFormatException {
        _reader = new CeosFileReader(leaderStream);
        _leaderFDR = new BaseRecord(_reader, -1, mission, leader_recordDefinitionFile);
        _reader.seek(_leaderFDR.getAbsolutPosition(_leaderFDR.getRecordLength()));
        _sceneHeaderRecord = new BaseSceneHeaderRecord(_reader, -1, mission, scene_recordDefinitionFile);
        _reader.seek(_sceneHeaderRecord.getAbsolutPosition(_sceneHeaderRecord.getRecordLength()));
        _mapProjRecord = new BaseRecord(_reader, -1, mission, mapproj_recordDefinitionFile);
        _reader.seek(_mapProjRecord.getAbsolutPosition(_mapProjRecord.getRecordLength()));
        _platformPositionRecord = new BaseRecord(_reader, -1, mission, platform_recordDefinitionFile);
        _reader.seek(_platformPositionRecord.getAbsolutPosition(_platformPositionRecord.getRecordLength()));
        _facilityRecord = new BaseRecord(_reader, -1, mission, facility_recordDefinitionFile);
    }

    public String getProductLevel() {
        return _sceneHeaderRecord.getAttributeString("Scene reference number");
    }

    public Calendar getDateImageWasTaken() {
        return _sceneHeaderRecord.getDateImageWasTaken();
    }

     /*
    public String getProcessingCode() throws IOException,
                                             IllegalCeosFormatException {
        return _sceneHeaderRecord.getProcessingCode();
    }

    public Calendar getDateImageWasTaken() {
        return _sceneHeaderRecord.getDateImageWasTaken();
    }    */

    public float[] getLatCorners() throws IOException,
                                           IllegalCeosFormatException {
        final double latUL = _mapProjRecord.getAttributeDouble("1st line 1st pixel geodetic latitude");
        final double latUR = _mapProjRecord.getAttributeDouble("1st line last valid pixel geodetic latitude");
        final double latLL = _mapProjRecord.getAttributeDouble("Last line 1st pixel geodetic latitude");
        final double latLR = _mapProjRecord.getAttributeDouble("Last line last valid pixel geodetic latitude");
        return new float[]{(float)latUL, (float)latUR, (float)latLL, (float)latLR};
    }

    public float[] getLonCorners() throws IOException,
                                           IllegalCeosFormatException {
        final double lonUL = _mapProjRecord.getAttributeDouble("1st line 1st pixel longitude");
        final double lonUR = _mapProjRecord.getAttributeDouble("1st line last valid pixel longitude");
        final double lonLL = _mapProjRecord.getAttributeDouble("Last line 1st pixel longitude");
        final double lonLR = _mapProjRecord.getAttributeDouble("Last line last valid pixel longitude");
        return new float[]{(float)lonUL, (float)lonUR, (float)lonLL, (float)lonLR};
    }

    public String getUsedProjection() throws IOException,
                                             IllegalCeosFormatException {
        return _mapProjRecord.getAttributeString("Map projection descriptor");
       // return _sceneHeaderRecord.getMapProjectionMethod().trim();
//        if (PROJECTION_KEY_RAW.equals(projKey)) {
//            return PalsarConstants.MAP_PROJECTION_RAW;
//        } else if (PROJECTION_KEY_UTM.equals(projKey)) {
//            return PalsarConstants.MAP_PROJECTION_UTM;
//        } else if (PROJECTION_KEY_PS.equals(projKey)) {
//            return PalsarConstants.MAP_PROJECTION_PS;
//        }
//        return PalsarConstants.MAP_PROJECTION_UNKNOWN;
    }

  /*  public long getNominalPixelsPerLine_1A_1B1() throws IOException,
                                                        IllegalCeosFormatException {
        return _ancillary1Record.getNumNominalPixelsPerLine_1A_1B1();
    }

    public long getNominalLinesPerScene_1A_1B1() throws IOException,
                                                        IllegalCeosFormatException {
        return _ancillary1Record.getNumNominalLinesPerScene_1A_1B1();
    }

    public double getNominalInterPixelDistance_1A_1B1() throws IOException,
                                                               IllegalCeosFormatException {
        return _ancillary1Record.getNominalInterPixelDistance_1A_1B1();
    }

    public double getNominalInterLineDistance_1A_1B1() throws IOException,
                                                              IllegalCeosFormatException {
        return _ancillary1Record.getNominalInterLineDistance_1A_1B1();
    }

    public double getImageSkew() throws IOException,
                                        IllegalCeosFormatException {
        return _ancillary1Record.getImageSkew();
    }

    private double getNumNominalPixelsPerLine() throws IOException,
                                                       IllegalCeosFormatException {
        return _ancillary1Record.getNumNominalPixelsPerLine();
    }

    private double getNumNominalLinesPerScene() throws IOException,
                                                       IllegalCeosFormatException {
        return _ancillary1Record.getNumNominalLinesPerScene();
    }     */

    /**
     * Gets the pixel size in x direction in meters.
     *
     * @return the pixel size in x direction
     *
     * @throws IOException
     * @throws IllegalCeosFormatException
     */
  /*  public double getNominalInterPixelDistance() throws IOException,
                                                        IllegalCeosFormatException {
        return _ancillary1Record.getNominalInterPixelDistance();
    }  */

    /**
     * Gets the pixel size in y direction in meters.
     *
     * @return the pixel size in y direction
     *
     * @throws IOException
     * @throws IllegalCeosFormatException
     */
  /*  public double getNominalInterLineDistance() throws IOException,
                                                       IllegalCeosFormatException {
        return _ancillary1Record.getNominalInterLineDistance();
    }

    public String getReferenceEllipsoidName() throws IOException,
                                                     IllegalCeosFormatException {
        return _ancillary1Record.getReferenceEllipsoid();
    }

    public double getSemiMinorAxis() throws IOException,
                                            IllegalCeosFormatException {
        return _ancillary1Record.getSemiminorAxisOfReferenceEllipsoid();
    }

    public double getSemiMajorAxis() throws IOException,
                                            IllegalCeosFormatException {
        return _ancillary1Record.getSemimajorAxisOfReferenceEllipsoid();
    }

    public String getDatumName() throws IOException,
                                        IllegalCeosFormatException {
        return _ancillary1Record.getGeodeticCoordinateName();
    }

    public long getUTMZoneIndex() throws IOException,
                                         IllegalCeosFormatException {
        return _ancillary1Record.getUTMZoneNumber();
    }

    public boolean isUTMSouthHemisphere() throws IOException,
                                                 IllegalCeosFormatException {
        return _ancillary1Record.getHemisphere() == 1;
    }    */

    /**
     * Gets the easting of the scene center in kilometers.
     *
     * @return the easting of the scene center
     *
     * @throws IOException
     * @throws IllegalCeosFormatException
     */
  /*  public double getUTMEasting() throws IOException,
                                         IllegalCeosFormatException {
        return _ancillary1Record.getSceneCenterEasting();
    }   */

    /**
     * Gets the northing  of the scene center in kilometers.
     *
     * @return the northing of the scene center
     *
     * @throws IOException
     * @throws IllegalCeosFormatException
     */
  /*  public double getUTMNorthing() throws IOException,
                                          IllegalCeosFormatException {
        return _ancillary1Record.getSceneCenterNorthing();
    }    */


    /**
     * Gets the orientation angle of the UTM projection in degree.
     *
     * @return the orientation angle in degree
     *
     * @throws IOException
     * @throws IllegalCeosFormatException
     */
  /*  public double getUTMOrientationAngle() throws IOException,
                                                  IllegalCeosFormatException {
        final double radian = _ancillary1Record.getAngleBetweenMapUTMVerticalAndTrueNorth();
        return Math.toDegrees(radian);
    }

    public GeoPos getPSProjectionOrigin() throws IOException,
                                                 IllegalCeosFormatException {
        final double lat = _ancillary1Record.getMapProjOriginLat();
        final double lon = _ancillary1Record.getMapProjOriginLon();
        return new GeoPos((float) lat, (float) lon);
    }

    public GeoPos getPSReferencePoint() throws IOException,
                                               IllegalCeosFormatException {
        final double lat = _ancillary1Record.getPSReferenceLat();
        final double lon = _ancillary1Record.getPSReferenceLon();
        return new GeoPos((float) lat, (float) lon);

    }

    public double getPSXCoordinate() throws IOException,
                                            IllegalCeosFormatException {
        return _ancillary1Record.getSceneCenterX();
    }

    public double getPSYCoordinate() throws IOException,
                                            IllegalCeosFormatException {
        return _ancillary1Record.getSceneCenterY();
    }

    public double getPSOrientationAngle() throws IOException,
                                                 IllegalCeosFormatException {
        return _ancillary1Record.getAngleBetweenMapPSVerticalAndTrueNorth();
    }      */

    /**
     * Gets the coefficiants for corrected L1B2 products.
     *
     * @return An arrays of arrays in the following order:
     *         <ul>
     *         <li>[0][0..9] latitude coefficients</li>
     *         <li>[1][0..9] longitude coefficients</li>
     *         <li>[2][0..9] x coefficients</li>
     *         <li>[3][0..9] y coefficients</li>
     *         </ul>
     *
     * @throws IOException
     * @throws IllegalCeosFormatException
     */
  /*  public double[][] getCorrectedTransformationCoeffs() throws IOException,
                                                                IllegalCeosFormatException {
        return _ancillary1Record.getTransformationCoeffsL1B2();
    }

    public double[][] getUncorrectedTransformationCoeffs(final int bandIndex) throws
                                                                              IOException,
                                                                              IllegalCeosFormatException {
        return _ancillary1Record.getTransformationCoefficientsFor(bandIndex);
    }

    public double getAbsoluteCalibrationGain(final int bandIndex) throws IOException,
                                                                         IllegalCeosFormatException {
        return _ancillary2Record.getAbsoluteCalibrationGain(bandIndex);
    }

    public double getAbsoluteCalibrationOffset(final int bandIndex) throws IOException,
                                                                           IllegalCeosFormatException {
        return _ancillary2Record.getAbsoluteCalibrationOffset(bandIndex);
    }

    public double[] getAbsoluteCalibrationGains() throws IOException,
                                                         IllegalCeosFormatException {
        final double[] gains = new double[4];
        for (int i = 0; i < gains.length; i++) {
            gains[i] = getAbsoluteCalibrationGain(i + 1);
        }
        return gains;
    }

    public double[] getAbsoluteCalibrationOffsets() throws IOException,
                                                           IllegalCeosFormatException {
        final double[] offsets = new double[4];
        for (int i = 0; i < offsets.length; i++) {
            offsets[i] = getAbsoluteCalibrationOffset(i + 1);
        }
        return offsets;
    }

    public String getSensorMode() throws IOException,
                                         IllegalCeosFormatException {
        return _ancillary2Record.getSensorOperationMode();
    }

    public int getLowerLimitStrength() throws IOException,
                                              IllegalCeosFormatException {
        return _ancillary2Record.getLowerLimitOfStrengthAfterCorrection();
    }

    public int getUpperLimitStrength() throws IOException,
                                              IllegalCeosFormatException {
        return _ancillary2Record.getLowerLimitOfStrengthAfterCorrection();
    }

    public double[] getExposureCoefficients() throws IOException,
                                                     IllegalCeosFormatException {
        final double[] exposureCoeffs = new double[4];
        for (int i = 0; i < exposureCoeffs.length; i++) {
            exposureCoeffs[i] = _ancillary2Record.getExposureCoefficient(i + 1);
        }
        return exposureCoeffs;
    }

    public char[] getSensorGains() throws IOException,
                                          IllegalCeosFormatException {
        final String gains = _ancillary2Record.getSensorGains();
        return gains.toCharArray();
    }

    public double[] getDetectorTemperatures() throws IOException,
                                                     IllegalCeosFormatException {
        final double[] temperatures = new double[4];
        for (int i = 0; i < temperatures.length; i++) {
            temperatures[i] = _ancillary2Record.getDetectorTemperature(i + 1);
        }
        return temperatures;
    }

    public double[] getDetectorAssemblyTemperatures() throws IOException,
                                                             IllegalCeosFormatException {

        final double[] temperatures = new double[4];
        for (int i = 0; i < temperatures.length; i++) {
            temperatures[i] = _ancillary2Record.getDetectorAssemblyTemperature(i + 1);
        }
        return temperatures;
    }

    public double getSignalProcessingUnitTemperature() throws IOException,
                                                              IllegalCeosFormatException {

        return _ancillary2Record.getSignalProcessingUnitTemperature();
    }

    public int getNumEffectiveDataPoints() throws IOException,
                                                  IllegalCeosFormatException {
        return _ancillary3Record.getNumDataPoints();
    }

    public int getYearOfFirstPoint() throws IOException,
                                            IllegalCeosFormatException {
        return _ancillary3Record.getFirstPointYear();
    }

    public int getMonthOfFirstPoint() throws IOException,
                                             IllegalCeosFormatException {
        return _ancillary3Record.getFirstPointMonth();
    }

    public int getDayOfFirstPoint() throws IOException,
                                           IllegalCeosFormatException {
        return _ancillary3Record.getFirstPointDay();
    }

    public int getTotalDaysOfFirstPoint() throws IOException,
                                                 IllegalCeosFormatException {
        return _ancillary3Record.getFirstPointTotalDays();
    }

    public double getTotalSecondsOfFirstPoint() throws IOException,
                                                       IllegalCeosFormatException {
        return _ancillary3Record.getFirstPointTotalSeconds();
    }

    public double getPointsInterval() throws IOException,
                                             IllegalCeosFormatException {
        return _ancillary3Record.getIntervalTimeBetweenPoints();
    }

    public String getPlatformtReferenceCoordinateSystem() throws IOException,
                                                                 IllegalCeosFormatException {
        return _ancillary3Record.getReferenceCoordinateSystem();
    }

    public double getFlightDirectionPositionalError() throws IOException,
                                                             IllegalCeosFormatException {
        return _ancillary3Record.getPositionalErrorFlightDirection();
    }

    public double getFlightDirectionPositionaVerticallError() throws IOException,
                                                                     IllegalCeosFormatException {
        return _ancillary3Record.getPositionalErrorFlightVerticalDirection();
    }

    public double getRadiusDirectionPositionalError() throws IOException,
                                                             IllegalCeosFormatException {
        return _ancillary3Record.getPositionalErrorRadiusDirection();
    }

    public double getFlightDirectionVelocityError() throws IOException,
                                                           IllegalCeosFormatException {
        return _ancillary3Record.getVelocityErrorFlightDirection();
    }

    public double getFlightDirectionVelocityVerticalError() throws IOException,
                                                                   IllegalCeosFormatException {
        return _ancillary3Record.getVelocityErrorFlightVerticalDirection();
    }

    public double getRadiusDirectionVelocityError() throws IOException,
                                                           IllegalCeosFormatException {
        return _ancillary3Record.getVelocityErrorRadiusDirection();
    }

    public Ancillary3Record.DataPoint[] getDataPoints() throws IOException,
                                                               IllegalCeosFormatException {
        return _ancillary3Record.getDataPoints();
    }

    public boolean isLeapSecondUsed() throws IOException,
                                             IllegalCeosFormatException {
        return _ancillary3Record.getFlagLeapSecond() == 1;
    }

    public double[] getF4Coefficients() throws IOException,
                                               IllegalCeosFormatException {
        return _ancillary1Record.getF4FunctionCoeffs_1B2();
    }     */

    public void addLeaderMetadata(MetadataElement sphElem) {
        MetadataElement metadata = new MetadataElement("Leader File Descriptor");
         _leaderFDR.assignMetadataTo(metadata);
        sphElem.addElement(metadata);

        metadata = new MetadataElement("Scene Parameters");
        _sceneHeaderRecord.assignMetadataTo(metadata);
        sphElem.addElement(metadata);

        metadata = new MetadataElement("Map Projection");
        _mapProjRecord.assignMetadataTo(metadata);
        sphElem.addElement(metadata);

        metadata = new MetadataElement("Platform Position");
        _platformPositionRecord.assignMetadataTo(metadata);
        sphElem.addElement(metadata);

        metadata = new MetadataElement("Facility");
        _facilityRecord.assignMetadataTo(metadata);
        sphElem.addElement(metadata);
    }

    public MetadataElement getMapProjectionMetadata() throws IOException,
                                                             IllegalCeosFormatException {
        final MetadataElement projMetadata = new MetadataElement("Map Projection");

        addGeneralProjectionMetadata(projMetadata);

        return projMetadata;
    }

    public MetadataElement getRadiometricMetadata() throws IOException,
                                                           IllegalCeosFormatException {
        final MetadataElement radioMetadata = new MetadataElement("Radiometric Calibration");
      /*  addAttribute(radioMetadata, "SENSOR_MODE", ProductData.createInstance(getSensorMode()));
        addAttribute(radioMetadata, "LOWER_LIMIT_STRENGTH",
                     ProductData.createInstance(new int[]{getLowerLimitStrength()}));
        addAttribute(radioMetadata, "UPPER_LIMIT_STRENGTH",
                     ProductData.createInstance(new int[]{getUpperLimitStrength()}));
        addAttribute(radioMetadata, "EXPOSURE_COEFFICIENT_BAND", ProductData.createInstance(getExposureCoefficients()));
        final char[] gains = getSensorGains();
        addAttribute(radioMetadata, "SENSOR_GAIN_BAND.1",
                     ProductData.createInstance(gains[0] == ' ' ? "" : "Gain " + gains[0]));
        addAttribute(radioMetadata, "SENSOR_GAIN_BAND.2",
                     ProductData.createInstance(gains[1] == ' ' ? "" : "Gain " + gains[1]));
        addAttribute(radioMetadata, "SENSOR_GAIN_BAND.3",
                     ProductData.createInstance(gains[2] == ' ' ? "" : "Gain " + gains[2]));
        addAttribute(radioMetadata, "SENSOR_GAIN_BAND.4",
                     ProductData.createInstance(gains[3] == ' ' ? "" : "Gain " + gains[3]));
        addAttribute(radioMetadata, "DETECTOR_TEMPERATURE_BAND", ProductData.createInstance(getDetectorTemperatures()),
                     UNIT_DEGREE);
        addAttribute(radioMetadata, "DETECTOR_ASSEMBLY_TEMPERATURE_BAND",
                     ProductData.createInstance(getDetectorAssemblyTemperatures()), UNIT_DEGREE);
        addAttribute(radioMetadata, "SIGNAL_PROCESSING_UNIT_TEMPERATURE",
                     ProductData.createInstance(new double[]{getSignalProcessingUnitTemperature()}), UNIT_DEGREE);
        final double[] absGains = getAbsoluteCalibrationGains();
        final double[] absOffsets = getAbsoluteCalibrationOffsets();
        addAttribute(radioMetadata, "ABSOLUTE_GAIN_BAND", ProductData.createInstance(absGains));
        addAttribute(radioMetadata, "ABSOLUTE_OFFSET_BAND", ProductData.createInstance(absOffsets));
                      */
        return radioMetadata;
    }

    public MetadataElement getPlatformMetadata() throws IOException,
                                                        IllegalCeosFormatException {
        final MetadataElement platformMeta = new MetadataElement("Platform Position Data");

     /*   addAttribute(platformMeta, "NUMBER_EFFECTIVE_DATA_POINTS",
                     ProductData.createInstance(new int[]{getNumEffectiveDataPoints()}));
        addAttribute(platformMeta, "YEAR_OF_FIRST_POINT", ProductData.createInstance(new int[]{getYearOfFirstPoint()}));
        addAttribute(platformMeta, "MONTH_OF_FIRST_POINT",
                     ProductData.createInstance(new int[]{getMonthOfFirstPoint()}));
        addAttribute(platformMeta, "DAY_OF_FIRST_POINT", ProductData.createInstance(new int[]{getDayOfFirstPoint()}));
        addAttribute(platformMeta, "TOTAL_DAYS_OF_FIRST_POINT",
                     ProductData.createInstance(new int[]{getTotalDaysOfFirstPoint()}));
        addAttribute(platformMeta, "TOTAL_SECONDS_OF_FIRST_POINT",
                     ProductData.createInstance(new double[]{getTotalSecondsOfFirstPoint()}));
        addAttribute(platformMeta, "POINTS_INTERVAL_TIME",
                     ProductData.createInstance(new double[]{getPointsInterval()}), UNIT_SECOND);

        addAttribute(platformMeta, "REFERENCE_COORDINATE_SYSTEM",
                     ProductData.createInstance(getPlatformtReferenceCoordinateSystem()));

        addAttribute(platformMeta, "POSITIONAL_ERROR_FLIGHT_DIRECTION",
                     ProductData.createInstance(new double[]{getFlightDirectionPositionalError()}), UNIT_METER);
        addAttribute(platformMeta, "POSITIONAL_ERROR_VERTICAL_FLIGHT_DIRECTION",
                     ProductData.createInstance(new double[]{getFlightDirectionPositionaVerticallError()}), UNIT_METER);
        addAttribute(platformMeta, "POSITIONAL_ERROR_RADIUS_DIRECTION",
                     ProductData.createInstance(new double[]{getRadiusDirectionPositionalError()}),
                     UNIT_METER_PER_SECOND);
        addAttribute(platformMeta, "VELOCITY_ERROR_FLIGHT_DIRECTION",
                     ProductData.createInstance(new double[]{getFlightDirectionVelocityError()}),
                     UNIT_METER_PER_SECOND);
        addAttribute(platformMeta, "VELOCITY_ERROR_VERTICAL_FLIGHT_DIRECTION",
                     ProductData.createInstance(new double[]{getFlightDirectionVelocityVerticalError()}),
                     UNIT_METER_PER_SECOND);
        addAttribute(platformMeta, "VELOCITY_ERROR_RADIUS_DIRECTION",
                     ProductData.createInstance(new double[]{getRadiusDirectionVelocityError()}),
                     UNIT_DEGREE_PER_SECOND);

        final Ancillary3Record.DataPoint[] dataPoints = getDataPoints();
        for (int i = 0; i < dataPoints.length; i++) {
            final int pIndex = i + 1;
            final Ancillary3Record.DataPoint dataPoint = dataPoints[i];
            addAttribute(platformMeta, "DATA_POINT_" + pIndex + "_POSITIONAL_VECTOR_X",
                         ProductData.createInstance(new double[]{dataPoint.getPositionalVectorDataPointX()}));
            addAttribute(platformMeta, "DATA_POINT_" + pIndex + "_POSITIONAL_VECTOR_Y",
                         ProductData.createInstance(new double[]{dataPoint.getPositionalVectorDataPointY()}));
            addAttribute(platformMeta, "DATA_POINT_" + pIndex + "_POSITIONAL_VECTOR_Z",
                         ProductData.createInstance(new double[]{dataPoint.getPositionalVectorDataPointZ()}));
            addAttribute(platformMeta, "DATA_POINT_" + pIndex + "_VELOCITY_VECTOR_X",
                         ProductData.createInstance(new double[]{dataPoint.getVelocityVectorDataPointX()}));
            addAttribute(platformMeta, "DATA_POINT_" + pIndex + "_VELOCITY_VECTOR_Y",
                         ProductData.createInstance(new double[]{dataPoint.getVelocityVectorDataPointY()}));
            addAttribute(platformMeta, "DATA_POINT_" + pIndex + "_VELOCITY_VECTOR_Z",
                         ProductData.createInstance(new double[]{dataPoint.getVelocityVectorDataPointZ()}));
        }

        addAttribute(platformMeta, "LEAP_SECOND", ProductData.createInstance(String.valueOf(isLeapSecondUsed())));
           */
        return platformMeta;
    }

    public void close() throws IOException {
        _reader.close();
        _reader = null;
    }

    private void addGeneralProjectionMetadata(final MetadataElement projMeta) throws
                                                                              IOException,
                                                                              IllegalCeosFormatException {
     /*   addAttribute(projMeta, "REFERENCE_ELLIPSOID", ProductData.createInstance(getReferenceEllipsoidName()));
        addAttribute(projMeta, "SEMI_MAJOR_AXIS", ProductData.createInstance(new double[]{getSemiMajorAxis()}),
                     UNIT_METER);
        addAttribute(projMeta, "SEMI_MINOR_AXIS", ProductData.createInstance(new double[]{getSemiMinorAxis()}),
                     UNIT_METER);

        addAttribute(projMeta, "GEODETIC_DATUM", ProductData.createInstance(getDatumName()));   */

        final float[] latCorners = getLatCorners();
        final float[] lonCorners = getLonCorners();

        addAttribute(projMeta, "SCENE_UPPER_LEFT_LATITUDE", ProductData.createInstance(new float[]{latCorners[0]}),
                     UNIT_DEGREE);
        addAttribute(projMeta, "SCENE_UPPER_LEFT_LONGITUDE", ProductData.createInstance(new float[]{lonCorners[0]}),
                     UNIT_DEGREE);
        addAttribute(projMeta, "SCENE_UPPER_RIGHT_LATITUDE", ProductData.createInstance(new float[]{latCorners[1]}),
                     UNIT_DEGREE);
        addAttribute(projMeta, "SCENE_UPPER_RIGHT_LONGITUDE", ProductData.createInstance(new float[]{lonCorners[1]}),
                     UNIT_DEGREE);
        addAttribute(projMeta, "SCENE_LOWER_LEFT_LATITUDE", ProductData.createInstance(new float[]{latCorners[2]}),
                     UNIT_DEGREE);
        addAttribute(projMeta, "SCENE_LOWER_LEFT_LONGITUDE", ProductData.createInstance(new float[]{lonCorners[2]}),
                     UNIT_DEGREE);
        addAttribute(projMeta, "SCENE_LOWER_RIGHT_LATITUDE", ProductData.createInstance(new float[]{latCorners[3]}),
                     UNIT_DEGREE);
        addAttribute(projMeta, "SCENE_LOWER_RIGHT_LONGITUDE", ProductData.createInstance(new float[]{lonCorners[3]}),
                     UNIT_DEGREE);
    }

    private void addGeneralCorrectedMetadata(final MetadataElement projMeta) throws
                                                                             IllegalCeosFormatException,
                                                                             IOException {
   /*     addAttribute(projMeta, "PIXELS_PER_LINE",
                     ProductData.createInstance(new double[]{getNumNominalPixelsPerLine()}));
        addAttribute(projMeta, "LINES_PER_SCENE",
                     ProductData.createInstance(new double[]{getNumNominalLinesPerScene()}));
        addAttribute(projMeta, "PIXEL_SIZE_X_CENTER",
                     ProductData.createInstance(new double[]{getNominalInterPixelDistance()}), UNIT_METER);
        addAttribute(projMeta, "PIXEL_SIZE_Y_CENTER",
                     ProductData.createInstance(new double[]{getNominalInterLineDistance()}), UNIT_METER);  */
    }

    private static MetadataAttribute createAttribute(final String name, final ProductData data) {
        return new MetadataAttribute(name.toUpperCase(), data, true);
    }

    private static MetadataAttribute addAttribute(final MetadataElement platformMetadata, final String name,
                                           final ProductData data) {
        return addAttribute(platformMetadata, name, data, null);
    }

    private static MetadataAttribute addAttribute(final MetadataElement platformMetadata, final String name,
                                           final ProductData data, final String unit) {
        final MetadataAttribute attribute = createAttribute(name, data);
        if (unit != null) {
            attribute.setUnit(unit);
        }

        platformMetadata.addAttribute(attribute);
        return attribute;
    }
}
