/*
 * $Id: PrismProductDirectory.java,v 1.1 2009-04-28 14:37:13 lveci Exp $
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
package org.esa.beam.dataio.ceos.prism;

import org.esa.beam.dataio.ceos.CeosHelper;
import org.esa.beam.dataio.ceos.IllegalCeosFormatException;
import org.esa.beam.dataio.ceos.prism.records.PrismAncillary2Record;
import org.esa.beam.dataio.ceos.prism.records.SceneHeaderRecord;
import org.esa.beam.dataio.ceos.records.Ancillary1Record;
import org.esa.beam.dataio.ceos.records.Ancillary3Record;
import org.esa.beam.framework.datamodel.*;
import org.esa.beam.framework.dataop.maptransf.*;
import org.esa.beam.util.Debug;
import org.esa.beam.util.Guardian;
import org.esa.beam.util.TreeNode;
import org.esa.beam.util.StringUtils;
import org.esa.beam.util.math.FXYSum;

import javax.imageio.stream.FileImageInputStream;
import javax.imageio.stream.ImageInputStream;
import java.io.File;
import java.io.IOException;
import java.util.Calendar;

class PrismProductDirectory {

    private static final int METER_PER_KILOMETER = 1000;

    private final static String MAP_PROJECTION_CODE_RAW = "NNNNN";
    private final static String MAP_PROJECTION_CODE_UTM = "YNNNN";
    private final static String MAP_PROJECTION_CODE_PS = "NNNNY";

    private final File _baseDir;
    private PrismVolumeDirectoryFile _volumeDirectoryFile;
    private PrismImageFile[] _imageFiles;
    private PrismLeaderFile _leaderFile;
    private PrismTrailerFile _trailerFile;
    private PrismSupplementalFile _supplementalFile;

    private static final String UNIT_METER = "meter";
    private static final String UNIT_KILOMETER = "kilometer";
    private static final String UNIT_DEGREE = "degree";
    private static final String UNIT_SECOND = "second";
    private static final String UNIT_METER_PER_SECOND = "m/sec";
    private static final String UNIT_DEGREE_PER_SECOND = "deg/sec";

    PrismProductDirectory(final File dir) throws IOException,
            IllegalCeosFormatException {
        Guardian.assertNotNull("dir", dir);

        _baseDir = dir;
        _volumeDirectoryFile = new PrismVolumeDirectoryFile(_baseDir);
        _leaderFile = new PrismLeaderFile(createInputStream(_volumeDirectoryFile.getLeaderFileName()));
        _trailerFile = new PrismTrailerFile(createInputStream(_volumeDirectoryFile.getTrailerFileName()));
        if (!_leaderFile.getProductLevel().equalsIgnoreCase(PrismConstants.PRODUCT_LEVEL_1B2)) {
            _supplementalFile = new PrismSupplementalFile(
                    createInputStream(_volumeDirectoryFile.getSupplementalFileName()));
        }

        final String[] imageFileNames = _volumeDirectoryFile.getImageFileNames();
        _imageFiles = new PrismImageFile[imageFileNames.length];
        for (int i = 0; i < imageFileNames.length; i++) {
            _imageFiles[i] = new PrismImageFile(createInputStream(imageFileNames[i]));
        }
    }

    Product createProduct() throws IOException,
            IllegalCeosFormatException {
        final String productName = _leaderFile.getProductName();
        final String productType = _leaderFile.getProductType();
        int width = 0;
        final int overlap = _imageFiles.length > 1 ? 32 : 0;
        for (int i = 0; i < _imageFiles.length; i++) {
            final PrismImageFile imageFile = _imageFiles[i];
            width += imageFile.getWidth() - overlap;
        }
        final int sceneWidth = width;
//        final int sceneWidth = _leaderFile.getSceneWidth();
        final int sceneHeight = _leaderFile.getSceneHeight();
        final Product product = new Product(productName, productType, sceneWidth, sceneHeight);

        product.setFileLocation(_baseDir);

        addBand(product);

        final Calendar imageStartDate = _leaderFile.getDateImageWasTaken();
        imageStartDate.add(Calendar.MILLISECOND, _imageFiles[0].getTotalMillisInDayOfLine(0));
        final ProductData.UTC utcScanStartTime = ProductData.UTC.create(imageStartDate.getTime(),
                _imageFiles[0].getMicrosecondsOfLine(0));

        final Calendar imageEndDate = _leaderFile.getDateImageWasTaken();
        imageEndDate.add(Calendar.MILLISECOND, _imageFiles[0].getTotalMillisInDayOfLine(sceneHeight - 1));
        final ProductData.UTC utcScanStopTime = ProductData.UTC.create(imageEndDate.getTime(),
                _imageFiles[0].getMicrosecondsOfLine(
                        sceneHeight - 1));

        product.setStartTime(utcScanStartTime);
        product.setEndTime(utcScanStopTime);
        product.setDescription("PRISM product Level " + productType);

        addGeoCoding(product);

        addMetadataTo(product);

        return product;
    }

    private void addBand(final Product product) {
        final String bandName = "radiance_1"; // constant

        final Band band = new Band(bandName, ProductData.TYPE_UINT8,
                product.getSceneRasterWidth(),
                product.getSceneRasterHeight());
        band.setSpectralBandIndex(0); // constant
        band.setSpectralWavelength(645); // constant
        band.setSpectralBandwidth(250); // constant
        band.setUnit("mw / (m^2*sr*nm)"); // constant
        band.setScalingFactor(_leaderFile.getAncillary2Record().getAbsoluteCalibrationGain());
        band.setScalingOffset(_leaderFile.getAncillary2Record().getAbsoluteCalibrationOffset());
        band.setNoDataValue(0.00);
        band.setNoDataValueUsed(true);

        product.addBand(band);
    }

    void close() throws IOException {
        for (int i = 0; i < _imageFiles.length; i++) {
            _imageFiles[i].close();
        }
        _leaderFile.close();
        if (_supplementalFile != null) {
            _supplementalFile.close();
        }
        _trailerFile.close();
        _volumeDirectoryFile.close();
    }

    PrismImageFile[] getImageFiles() {
        return _imageFiles;
    }


    TreeNode<File> getProductComponents() {
        try {
            File parentDir = _baseDir.getParentFile();
            final TreeNode<File> root = new TreeNode<File>(parentDir.getCanonicalPath());
            root.setContent(parentDir);

            final TreeNode<File> dir = new TreeNode<File>(_baseDir.getName());
            dir.setContent(_baseDir);
            root.addChild(dir);

            File volumeFile = CeosHelper.getVolumeFile(_baseDir);
            addFileToDir(dir, volumeFile);

            addFileToDir(dir, new File(_volumeDirectoryFile.getLeaderFileName()));
            addFileToDir(dir, new File(_volumeDirectoryFile.getTrailerFileName()));
            String supplemental = _volumeDirectoryFile.getSupplementalFileName();
            if (StringUtils.isNotNullAndNotEmpty(supplemental)) {
                addFileToDir(dir, new File(supplemental));
            }
            final String[] imageFileNames = _volumeDirectoryFile.getImageFileNames();
            for (int i = 0; i < imageFileNames.length; i++) {
                addFileToDir(dir, new File(imageFileNames[i]));
            }

            return root;
        } catch (IOException e) {
            return null;
        }
    }

    private void addFileToDir(TreeNode<File> dir, File file) {
        final TreeNode<File> fileNode = new TreeNode<File>(file.getName());
        fileNode.setContent(file);
        dir.addChild(fileNode);
    }

    private void addGeoCoding(final Product product) throws IllegalCeosFormatException,
            IOException {
        final String projectionCode = getProjectionCode();
        if (MAP_PROJECTION_CODE_RAW.equalsIgnoreCase(projectionCode)) {
            final int overlap = 32;
            final int pixelOffsetX = overlap / 2;
            final PrismImageFile[] imageFiles = getImageFiles();

            final CombinedFXYGeoCoding.CodingWrapper[] codingWrappers;
            codingWrappers = new CombinedFXYGeoCoding.CodingWrapper[imageFiles.length];

            for (int i = 0; i < imageFiles.length; i++) {
                final PrismImageFile imageFile = imageFiles[i];
                int width = imageFile.getWidth() - overlap;
                final int height = imageFile.getHeight();
                final int offset = width * i;
                final int ccdNumber = imageFile.getImageNumber();

                final double[][] uncorrectedCoeffs = _leaderFile.getUncorrectedTransformationCoeffs(ccdNumber);


                final FXYSum.Cubic funcLat = new FXYSum.Cubic(CeosHelper.sortToFXYSumOrder(uncorrectedCoeffs[0]));
                final FXYSum.Cubic funcLon = new FXYSum.Cubic(CeosHelper.sortToFXYSumOrder(uncorrectedCoeffs[1]));
                final FXYSum.Cubic funcX = new FXYSum.Cubic(CeosHelper.sortToFXYSumOrder(uncorrectedCoeffs[2]));
                final FXYSum.Cubic funcY = new FXYSum.Cubic(CeosHelper.sortToFXYSumOrder(uncorrectedCoeffs[3]));

                final FXYGeoCoding gc = new FXYGeoCoding(pixelOffsetX, 0.0f, 1.0f, 1.0f,
                        funcX, funcY,
                        funcLat, funcLon,
                        Datum.ITRF_97);
                if (i != imageFiles.length - 1) {
                    ++width;
                }
                codingWrappers[i] = new CombinedFXYGeoCoding.CodingWrapper(gc, offset, 0, width, height);
            }
            final GeoCoding gc = new CombinedFXYGeoCoding(codingWrappers);
            product.setGeoCoding(gc);

        } else if (MAP_PROJECTION_CODE_UTM.equalsIgnoreCase(projectionCode)) {
            final int meterPerKilometer = 1000;
            final int sceneWidth = product.getSceneRasterWidth();
            final int sceneHeight = product.getSceneRasterHeight();

            final int zoneIndex = (int) _leaderFile.getUTMZoneIndex();
            final boolean isSouth = _leaderFile.isUTMSouthHemisphere();
            double easting = _leaderFile.getUTMEasting() * meterPerKilometer; // km -> meter
            double northing = _leaderFile.getUTMNorthing() * meterPerKilometer;    // km -> meter

            final double pixelSizeX = _leaderFile.getNominalInterPixelDistance();
            final double pixelSizeY = _leaderFile.getNominalInterLineDistance();
            final float orientationAngle = (float) _leaderFile.getUTMOrientationAngle();

            final MapInfo mapInfo = new MapInfo(UTM.createProjection(zoneIndex - 1, isSouth),
                    sceneWidth * 0.5f, sceneHeight * 0.5f,
                    (float) easting, (float) northing,
                    (float) pixelSizeX, (float) pixelSizeY, Datum.ITRF_97);
            // the BEAM convention for rotation angle uses opposite sign (rq - 16.10.2008)
            mapInfo.setOrientation(-orientationAngle);
            mapInfo.setSceneWidth(sceneWidth);
            mapInfo.setSceneHeight(sceneHeight);
            product.setGeoCoding(new MapGeoCoding(mapInfo));

        } else if (MAP_PROJECTION_CODE_PS.equalsIgnoreCase(projectionCode)) {
            final double[] parameterValues = StereographicDescriptor.PARAMETER_DEFAULT_VALUES;
            parameterValues[0] = Ellipsoid.GRS_80.getSemiMajor();
            parameterValues[1] = Ellipsoid.GRS_80.getSemiMinor();
            final GeoPos psReferencePoint = _leaderFile.getPSReferencePoint();
            final GeoPos psProjectionOrigin = _leaderFile.getPSProjectionOrigin();

            parameterValues[2] = psProjectionOrigin.getLat();         // Latitude_True_Scale
            parameterValues[3] = psReferencePoint.getLon();       // Central_Meridian


            final MapTransform transform = MapTransformFactory.createTransform(StereographicDescriptor.TYPE_ID,
                    parameterValues);
            final MapProjection projection = new MapProjection(StereographicDescriptor.NAME, transform);
            final double pixelSizeX = _leaderFile.getNominalInterPixelDistance();
            final double pixelSizeY = _leaderFile.getNominalInterLineDistance();
            final double easting = _leaderFile.getPSXCoordinate() * METER_PER_KILOMETER;
            final double northing = _leaderFile.getPSYCoordinate() * METER_PER_KILOMETER;
            final int sceneRasterWidth = product.getSceneRasterWidth();
            final int sceneRasterHeight = product.getSceneRasterHeight();
            final MapInfo mapInfo = new MapInfo(projection,
                    sceneRasterWidth * 0.5f, sceneRasterHeight * 0.5f,
                    (float) easting, (float) northing,
                    (float) pixelSizeX, (float) pixelSizeY, Datum.ITRF_97);
            mapInfo.setOrientation((float) _leaderFile.getPSOrientationAngle());
            mapInfo.setSceneWidth(sceneRasterWidth);
            mapInfo.setSceneHeight(sceneRasterHeight);
            product.setGeoCoding(new MapGeoCoding(mapInfo));

            // Alternative geo-coding for polar-stereographic
//            final double[][] l1B2Coeffs = _leaderFile.getCorrectedTransformationCoeffs();
//            final FXYSum.Cubic funcLat = new FXYSum.Cubic(CeosHelper.sortToFXYSumOrder(l1B2Coeffs[0]));
//            final FXYSum.Cubic funcLon = new FXYSum.Cubic(CeosHelper.sortToFXYSumOrder(l1B2Coeffs[1]));
//            final FXYSum.Cubic funcX = new FXYSum.Cubic(CeosHelper.sortToFXYSumOrder(l1B2Coeffs[2]));
//            final FXYSum.Cubic funcY = new FXYSum.Cubic(CeosHelper.sortToFXYSumOrder(l1B2Coeffs[3]));
//            final FXYGeoCoding funcGeoCoding = new FXYGeoCoding(0.0f, 0.0f, 1.0f, 1.0f,
//                                                                funcX, funcY, funcLat, funcLon,
//                                                                Datum.ITRF_97);
//            product.setGeoCoding(funcGeoCoding);

//            throw new IOException("Polar stereografic projection not implemented for PRISM products");
        } else {
            Debug.trace("Unknown map projection method. Could not create geo-coding.");
        }
    }

    private void addMetadataTo(final Product product)
            throws
            IOException,
            IllegalCeosFormatException {
        final MetadataElement metadataRoot = product.getMetadataRoot();
        assignSPH(metadataRoot);
        _volumeDirectoryFile.assignMetadataTo(metadataRoot);
    }

    private void assignSPH(MetadataElement metadataRoot) throws IOException, IllegalCeosFormatException {
        final MetadataElement sphElement = new MetadataElement("SPH");
        metadataRoot.addElement(sphElement);
        sphElement.addElement(getMapProjectionMetadata());
        sphElement.addElement(getRadiometricMetadata());
        sphElement.addElement(getPlatformMetadata());
    }

    private MetadataElement getMapProjectionMetadata() throws IOException,
            IllegalCeosFormatException {
        final MetadataElement projMetadata = new MetadataElement("Map Projection");

        addGeneralProjectionMetadata(projMetadata);

        final String usedProjection = getProjectionCode();
        if (usedProjection.equalsIgnoreCase(PrismProductDirectory.MAP_PROJECTION_CODE_RAW)) {
            addRawProjectionMetadata(projMetadata);
        } else if (usedProjection.equalsIgnoreCase(PrismProductDirectory.MAP_PROJECTION_CODE_UTM)) {
            addGeneralCorrectedMetadata(projMetadata);
            addUTMProjectionMetadata(projMetadata);
        } else if (usedProjection.equalsIgnoreCase(PrismProductDirectory.MAP_PROJECTION_CODE_PS)) {
            addGeneralCorrectedMetadata(projMetadata);
            addPSProjectionMetadata(projMetadata);
        }
        return projMetadata;
    }


    private void addGeneralProjectionMetadata(final MetadataElement projMeta) throws
            IOException,
            IllegalCeosFormatException {
        final PrismLeaderFile lf = _leaderFile;
        final Ancillary1Record a1r = lf.getAncillary1Record();
        addAttribute(projMeta, "REFERENCE_ELLIPSOID", ProductData.createInstance(a1r.getReferenceEllipsoid()));
        addAttribute(projMeta, "SEMI_MAJOR_AXIS", ProductData.createInstance(new double[]{lf.getSemiMajorAxis()}),
                UNIT_METER);
        addAttribute(projMeta, "SEMI_MINOR_AXIS", ProductData.createInstance(new double[]{lf.getSemiMinorAxis()}),
                UNIT_METER);

        addAttribute(projMeta, "GEODETIC_DATUM", ProductData.createInstance(lf.getDatumName()));

        final double[] latCorners = getLatCorners();
        final double[] lonCorners = getLonCorners();

        addAttribute(projMeta, "SCENE_UPPER_LEFT_LATITUDE", ProductData.createInstance(new double[]{latCorners[0]}),
                UNIT_DEGREE);
        addAttribute(projMeta, "SCENE_UPPER_LEFT_LONGITUDE", ProductData.createInstance(new double[]{lonCorners[0]}),
                UNIT_DEGREE);
        addAttribute(projMeta, "SCENE_UPPER_RIGHT_LATITUDE", ProductData.createInstance(new double[]{latCorners[1]}),
                UNIT_DEGREE);
        addAttribute(projMeta, "SCENE_UPPER_RIGHT_LONGITUDE", ProductData.createInstance(new double[]{lonCorners[1]}),
                UNIT_DEGREE);
        addAttribute(projMeta, "SCENE_LOWER_LEFT_LATITUDE", ProductData.createInstance(new double[]{latCorners[2]}),
                UNIT_DEGREE);
        addAttribute(projMeta, "SCENE_LOWER_LEFT_LONGITUDE", ProductData.createInstance(new double[]{lonCorners[2]}),
                UNIT_DEGREE);
        addAttribute(projMeta, "SCENE_LOWER_RIGHT_LATITUDE", ProductData.createInstance(new double[]{latCorners[3]}),
                UNIT_DEGREE);
        addAttribute(projMeta, "SCENE_LOWER_RIGHT_LONGITUDE", ProductData.createInstance(new double[]{lonCorners[3]}),
                UNIT_DEGREE);
    }

    private void addRawProjectionMetadata(final MetadataElement projMeta) throws
            IOException,
            IllegalCeosFormatException {
        final PrismLeaderFile lf = _leaderFile;
        for (int i = 1; i <= 4; i++) {
            final double[][] uncorrectedTransformationCoeffs = lf.getUncorrectedTransformationCoeffs(i);

            for (int j = 0; j < uncorrectedTransformationCoeffs[0].length; j++) {
                final double coeffLat = uncorrectedTransformationCoeffs[0][j];
                addAttribute(projMeta, "CCD[" + i + "]_COEFFICIENTS_LATITUDE." + j,
                        ProductData.createInstance(new double[]{coeffLat}));
            }
            for (int j = 0; j < uncorrectedTransformationCoeffs[1].length; j++) {
                final double coeffLon = uncorrectedTransformationCoeffs[1][j];
                addAttribute(projMeta, "CCD[" + i + "]_COEFFICIENTS_LONGITUDE." + j,
                        ProductData.createInstance(new double[]{coeffLon}));
            }
            for (int j = 0; j < uncorrectedTransformationCoeffs[2].length; j++) {
                final double coeffX = uncorrectedTransformationCoeffs[2][j];
                addAttribute(projMeta, "CCD[" + i + "]_COEFFICIENTS_X." + j,
                        ProductData.createInstance(new double[]{coeffX}));
            }
            for (int j = 0; j < uncorrectedTransformationCoeffs[3].length; j++) {
                final double coeffY = uncorrectedTransformationCoeffs[3][j];
                addAttribute(projMeta, "CCD[" + i + "]_COEFFICIENTS_Y." + j,
                        ProductData.createInstance(new double[]{coeffY}));
            }
        }

        addAttribute(projMeta, "PIXELS_PER_LINE",
                ProductData.createInstance(new long[]{lf.getNominalPixelsPerLine_1A_1B1()}));
        addAttribute(projMeta, "LINES_PER_SCENE",
                ProductData.createInstance(new long[]{lf.getNominalLinesPerScene_1A_1B1()}));
        addAttribute(projMeta, "PIXEL_SIZE_X_CENTER",
                ProductData.createInstance(new double[]{lf.getNominalInterPixelDistance_1A_1B1()}), UNIT_METER);
        addAttribute(projMeta, "PIXEL_SIZE_Y_CENTER",
                ProductData.createInstance(new double[]{lf.getNominalInterLineDistance_1A_1B1()}), UNIT_METER);
        addAttribute(projMeta, "IMAGE_SKEW_CENTER", ProductData.createInstance(new double[]{lf.getImageSkew()}),
                "milliradian");

    }

    private void addUTMProjectionMetadata(final MetadataElement projMeta) throws
            IOException,
            IllegalCeosFormatException {
        final PrismLeaderFile lf = _leaderFile;
        addAttribute(projMeta, "HEMISPHERE", ProductData.createInstance(lf.isUTMSouthHemisphere() ? "South" : "North"));
        addAttribute(projMeta, "UTM_ZONE_NUMBER", ProductData.createInstance(new long[]{lf.getUTMZoneIndex()}));
        addAttribute(projMeta, "UTM_NORTHING", ProductData.createInstance(new double[]{lf.getUTMNorthing()}),
                UNIT_KILOMETER);
        addAttribute(projMeta, "UTM_EASTING", ProductData.createInstance(new double[]{lf.getUTMEasting()}),
                UNIT_KILOMETER);
        final MetadataAttribute orientation = addAttribute(projMeta, "ORIENTATION",
                ProductData.createInstance(new double[]{
                        lf.getUTMOrientationAngle()
                }),
                UNIT_DEGREE);
        orientation.setDescription("Angle between the map projection vertical axis and the true north at scene center");

    }

    private void addPSProjectionMetadata(final MetadataElement projMeta) throws
            IOException,
            IllegalCeosFormatException {
        final PrismLeaderFile lf = _leaderFile;
        final GeoPos origin = lf.getPSProjectionOrigin();
        addAttribute(projMeta, "MAP_PROJECTION_ORIGIN",
                ProductData.createInstance(origin.getLatString() + " , " + origin.getLonString()));

        final GeoPos reference = lf.getPSReferencePoint();
        addAttribute(projMeta, "REFERENCE_POINT",
                ProductData.createInstance(reference.getLatString() + " , " + reference.getLonString()));

        addAttribute(projMeta, "COORDINATE_CENTER_X", ProductData.createInstance(new double[]{lf.getPSXCoordinate()}),
                UNIT_KILOMETER);
        addAttribute(projMeta, "COORDINATE_CENTER_Y)", ProductData.createInstance(new double[]{lf.getPSYCoordinate()}),
                UNIT_KILOMETER);

        final MetadataAttribute orientation = addAttribute(projMeta, "ORIENTATION",
                ProductData.createInstance(
                        new double[]{lf.getPSOrientationAngle()}),
                UNIT_DEGREE);
        orientation.setDescription("Angle between the map projection vertical axis and the true north at scene center");
    }

    private void addGeneralCorrectedMetadata(final MetadataElement projMeta) throws
            IllegalCeosFormatException,
            IOException {
        final Ancillary1Record ar = _leaderFile.getAncillary1Record();
        addAttribute(projMeta, "PIXELS_PER_LINE", ProductData.createInstance(
                new double[]{ar.getNumNominalPixelsPerLine()}));
        addAttribute(projMeta, "LINES_PER_SCENE", ProductData.createInstance(
                new double[]{ar.getNumNominalLinesPerScene()}));
        addAttribute(projMeta, "PIXEL_SIZE_X_CENTER", ProductData.createInstance(
                new double[]{ar.getNominalInterPixelDistance()}), UNIT_METER);
        addAttribute(projMeta, "PIXEL_SIZE_Y_CENTER", ProductData.createInstance(
                new double[]{ar.getNominalInterLineDistance()}), UNIT_METER);
    }

    private MetadataElement getRadiometricMetadata() {
        final PrismAncillary2Record ar = _leaderFile.getAncillary2Record();

        final MetadataElement radioMetadata = new MetadataElement("Radiometric Calibration");
        addAttribute(radioMetadata, "SENSOR_MODE", ProductData.createInstance(ar.getSensorOperationMode()));
        addAttribute(radioMetadata, "LOWER_LIMIT_STRENGTH",
                ProductData.createInstance(new int[]{ar.getLowerLimitOfStrengthAfterCorrection()}));
        addAttribute(radioMetadata, "UPPER_LIMIT_STRENGTH",
                ProductData.createInstance(new int[]{ar.getLowerLimitOfStrengthAfterCorrection()}));
        final char[] gains = ar.getSensorGains().trim().toCharArray();
        for (int i = 0; i < gains.length; i++) {
            addAttribute(radioMetadata, "SENSOR_GAIN." + (i + 1),
                    ProductData.createInstance(gains[i] == ' ' ? "" : "Gain " + gains[i]));
        }

        addAttribute(radioMetadata, "SIGNAL_PROCESSING_UNIT_TEMPERATURE",
                ProductData.createInstance(new double[]{ar.getSignalProcessingSectionTemperature()}), UNIT_DEGREE);
        final double absGain = ar.getAbsoluteCalibrationGain();
        final double absOffset = ar.getAbsoluteCalibrationOffset();
        addAttribute(radioMetadata, "ABSOLUTE_GAIN_BAND", ProductData.createInstance(new double[]{absGain}));
        addAttribute(radioMetadata, "ABSOLUTE_OFFSET_BAND", ProductData.createInstance(new double[]{absOffset}));

        return radioMetadata;
    }

    private MetadataElement getPlatformMetadata() {
        final Ancillary3Record ar = _leaderFile.getAncillary3Record();

        final MetadataElement platformMeta = new MetadataElement("Platform Position Data");

        addAttribute(platformMeta, "NUMBER_EFFECTIVE_DATA_POINTS",
                ProductData.createInstance(new int[]{ar.getNumDataPoints()}));
        addAttribute(platformMeta, "YEAR_OF_FIRST_POINT",
                ProductData.createInstance(new int[]{ar.getFirstPointYear()}));
        addAttribute(platformMeta, "MONTH_OF_FIRST_POINT",
                ProductData.createInstance(new int[]{ar.getFirstPointMonth()}));
        addAttribute(platformMeta, "DAY_OF_FIRST_POINT", ProductData.createInstance(new int[]{ar.getFirstPointDay()}));
        addAttribute(platformMeta, "TOTAL_DAYS_OF_FIRST_POINT",
                ProductData.createInstance(new int[]{ar.getFirstPointTotalDays()}));
        addAttribute(platformMeta, "TOTAL_SECONDS_OF_FIRST_POINT",
                ProductData.createInstance(new double[]{ar.getFirstPointTotalSeconds()}));
        addAttribute(platformMeta, "POINTS_INTERVAL_TIME",
                ProductData.createInstance(new double[]{ar.getIntervalTimeBetweenPoints()}), UNIT_SECOND);

        addAttribute(platformMeta, "REFERENCE_COORDINATE_SYSTEM",
                ProductData.createInstance(ar.getReferenceCoordinateSystem()));

        addAttribute(platformMeta, "POSITIONAL_ERROR_FLIGHT_DIRECTION",
                ProductData.createInstance(new double[]{ar.getPositionalErrorFlightDirection()}), UNIT_METER);
        addAttribute(platformMeta, "POSITIONAL_ERROR_VERTICAL_FLIGHT_DIRECTION",
                ProductData.createInstance(new double[]{ar.getPositionalErrorFlightVerticalDirection()}),
                UNIT_METER);
        addAttribute(platformMeta, "POSITIONAL_ERROR_RADIUS_DIRECTION",
                ProductData.createInstance(new double[]{ar.getPositionalErrorRadiusDirection()}),
                UNIT_METER_PER_SECOND);
        addAttribute(platformMeta, "VELOCITY_ERROR_FLIGHT_DIRECTION",
                ProductData.createInstance(new double[]{ar.getVelocityErrorFlightDirection()}),
                UNIT_METER_PER_SECOND);
        addAttribute(platformMeta, "VELOCITY_ERROR_VERTICAL_FLIGHT_DIRECTION",
                ProductData.createInstance(new double[]{ar.getVelocityErrorFlightVerticalDirection()}),
                UNIT_METER_PER_SECOND);
        addAttribute(platformMeta, "VELOCITY_ERROR_RADIUS_DIRECTION",
                ProductData.createInstance(new double[]{ar.getVelocityErrorRadiusDirection()}),
                UNIT_DEGREE_PER_SECOND);

        final Ancillary3Record.DataPoint[] dataPoints = ar.getDataPoints();
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

        addAttribute(platformMeta, "LEAP_SECOND",
                ProductData.createInstance(String.valueOf(ar.getFlagLeapSecond() == 1)));

        return platformMeta;
    }

    private double[] getLatCorners() {
        final SceneHeaderRecord shr = _leaderFile.getSceneHeaderRecord();
        final double latUL = shr.getSceneCornerUpperLeftLat();
        final double latUR = shr.getSceneCornerUpperRightLat();
        final double latLL = shr.getSceneCornerLowerLeftLat();
        final double latLR = shr.getSceneCornerLowerRightLat();
        return new double[]{latUL, latUR, latLL, latLR};
    }

    private double[] getLonCorners() {
        final SceneHeaderRecord shr = _leaderFile.getSceneHeaderRecord();
        final double lonUL = shr.getSceneCornerUpperLeftLon();
        final double lonUR = shr.getSceneCornerUpperRightLon();
        final double lonLL = shr.getSceneCornerLowerLeftLon();
        final double lonLR = shr.getSceneCornerLowerLeftLat();
        return new double[]{lonUL, lonUR, lonLL, lonLR};
    }

    private String getProjectionCode() {
        return _leaderFile.getSceneHeaderRecord().getMapProjectionMethod().trim();
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


    private ImageInputStream createInputStream(final String fileName) throws IOException {
        return new FileImageInputStream(new File(_baseDir, fileName));
    }
}
