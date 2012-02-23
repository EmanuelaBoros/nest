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
package org.esa.nest.dataio.ceos;

import org.esa.beam.framework.datamodel.MetadataElement;
import org.esa.nest.dataio.binary.BinaryRecord;
import org.esa.nest.dataio.binary.BinaryDBReader;
import org.esa.nest.dataio.binary.BinaryFileReader;

import javax.imageio.stream.ImageInputStream;
import java.io.IOException;


public class CEOSLeaderFile {

    protected BinaryRecord _leaderFDR = null;
    protected BinaryRecord _sceneHeaderRecord = null;
    protected BinaryRecord _platformPositionRecord = null;
    protected BinaryRecord _mapProjRecord = null;
    protected BinaryRecord _dataQualityRecord = null;
    protected BinaryRecord _histogramRecord = null;
    protected BinaryRecord _attitudeRecord = null;
    protected BinaryRecord _radiometricRecord = null;
    protected BinaryRecord _radiometricCompRecord = null;
    protected BinaryRecord _detailedProcessingRecord = null;
    protected BinaryRecord _facilityRecord = null;

    protected final static String scene_recordDefinitionFile = "scene_record.xml";
    protected final static String platformPosition_recordDefinitionFile = "platform_position_record.xml";
    protected final static String mapproj_recordDefinitionFile = "map_proj_record.xml";
    protected final static String dataQuality_recordDefinitionFile = "data_quality_summary_record.xml";
    protected final static String histogram_recordDefinitionFile = "data_histogram_record.xml";
    protected final static String attitude_recordDefinitionFile = "attitude_record.xml";
    protected final static String radiometric_recordDefinitionFile = "radiometric_record.xml";
    protected final static String radiometric_comp_recordDefinitionFile = "radiometric_compensation_record.xml";
    protected final static String detailedProcessing_recordDefinitionFile = "detailed_processing_record.xml";
    protected final static String facility_recordDefinitionFile = "facility_record.xml";

    private static org.jdom.Document sceneXML;
    private static org.jdom.Document mapProjXML;
    private static org.jdom.Document platformXML;
    private static org.jdom.Document attitudeXML;
    private static org.jdom.Document radiometricXML;
    private static org.jdom.Document radiometricCompXML;
    private static org.jdom.Document dataQualityXML;
    private static org.jdom.Document histogramXML;
    private static org.jdom.Document detailProcXML;
    private static org.jdom.Document facilityXML;

    protected CEOSLeaderFile() {
    }

    protected CEOSLeaderFile(final ImageInputStream stream, final String mission, final String defnFile)
            throws IOException {
        final BinaryFileReader reader = new BinaryFileReader(stream);

        org.jdom.Document leaderXML = BinaryDBReader.loadDefinitionFile(mission, defnFile);
        _leaderFDR = new BinaryRecord(reader, -1, leaderXML);
        reader.seek(_leaderFDR.getRecordEndPosition());

        if(sceneXML == null)
            sceneXML = BinaryDBReader.loadDefinitionFile(mission, scene_recordDefinitionFile);
        int num = _leaderFDR.getAttributeInt("Number of data set summary records");
        for(int i=0; i < num; ++i) {
            _sceneHeaderRecord = new BinaryRecord(reader, -1, sceneXML);
            reader.seek(_sceneHeaderRecord.getRecordEndPosition());
        }

        if(mapProjXML == null)
            mapProjXML = BinaryDBReader.loadDefinitionFile(mission, mapproj_recordDefinitionFile);
        num = _leaderFDR.getAttributeInt("Number of map projection data records");
        for(int i=0; i < num; ++i) {
            _mapProjRecord = new BinaryRecord(reader, -1, mapProjXML);
            reader.seek(_mapProjRecord.getRecordEndPosition());
        }

        if(platformXML == null)
            platformXML = BinaryDBReader.loadDefinitionFile(mission, platformPosition_recordDefinitionFile);
        num = _leaderFDR.getAttributeInt("Number of platform pos. data records");
        for(int i=0; i < num; ++i) {
            _platformPositionRecord = new BinaryRecord(reader, -1, platformXML);
            reader.seek(_platformPositionRecord.getRecordEndPosition());
        }

        if(attitudeXML == null)
            attitudeXML = BinaryDBReader.loadDefinitionFile(mission, attitude_recordDefinitionFile);
        num = _leaderFDR.getAttributeInt("Number of attitude data records");
        for(int i=0; i < num; ++i) {
            _attitudeRecord = new BinaryRecord(reader, -1, attitudeXML);
            reader.seek(_attitudeRecord.getRecordEndPosition());
        }

        if(radiometricXML == null)
            radiometricXML = BinaryDBReader.loadDefinitionFile(mission, radiometric_recordDefinitionFile);
        num = _leaderFDR.getAttributeInt("Number of radiometric data records");
        for(int i=0; i < num; ++i) {
            _radiometricRecord = new BinaryRecord(reader, -1, radiometricXML);
            reader.seek(_radiometricRecord.getRecordEndPosition());
        }

        if(radiometricCompXML == null)
            radiometricCompXML = BinaryDBReader.loadDefinitionFile(mission, radiometric_comp_recordDefinitionFile);
        num = _leaderFDR.getAttributeInt("Number of rad. compensation records");
        for(int i=0; i < num; ++i) {
            _radiometricCompRecord = new BinaryRecord(reader, -1, radiometricCompXML);
            reader.seek(_radiometricCompRecord.getRecordEndPosition());
        }

        if(dataQualityXML == null)
            dataQualityXML = BinaryDBReader.loadDefinitionFile(mission, dataQuality_recordDefinitionFile);
        num = _leaderFDR.getAttributeInt("Number of data quality summary records");
        for(int i=0; i < num; ++i) {
            _dataQualityRecord = new BinaryRecord(reader, -1, dataQualityXML);
            reader.seek(_dataQualityRecord.getRecordEndPosition());
        }

        if(histogramXML == null)
            histogramXML = BinaryDBReader.loadDefinitionFile(mission, histogram_recordDefinitionFile);
        num = _leaderFDR.getAttributeInt("Number of data histograms records");
        for(int i=0; i < num; ++i) {
            _histogramRecord = new BinaryRecord(reader, -1, histogramXML);
            reader.seek(_histogramRecord.getRecordEndPosition());
        }

        if(detailProcXML == null)
            detailProcXML = BinaryDBReader.loadDefinitionFile(mission, detailedProcessing_recordDefinitionFile);
        num = _leaderFDR.getAttributeInt("Number of det. processing records");
        for(int i=0; i < num; ++i) {
            _detailedProcessingRecord = new BinaryRecord(reader, -1, detailProcXML);
            reader.seek(_detailedProcessingRecord.getRecordEndPosition());
        }

        if(facilityXML == null)
            facilityXML = BinaryDBReader.loadDefinitionFile(mission, facility_recordDefinitionFile);
        num = _leaderFDR.getAttributeInt("Number of facility data records");
        for(int i=0; i < num; ++i) {
            _facilityRecord = new BinaryRecord(reader, -1, facilityXML);
            reader.seek(_facilityRecord.getRecordEndPosition());
        }

        reader.close();
    }

    public final BinaryRecord getSceneRecord() {
        return _sceneHeaderRecord;
    }

    public BinaryRecord getMapProjRecord() {
        return _mapProjRecord;
    }

    public BinaryRecord getPlatformPositionRecord() {
        return _platformPositionRecord;
    }

    public BinaryRecord getHistogramRecord() {
        return _histogramRecord;
    }

    public BinaryRecord getRadiometricRecord() {
        return _radiometricRecord;
    }

    public BinaryRecord getFacilityRecord() {
        return _facilityRecord;
    }

    public BinaryRecord getDetailedProcessingRecord() {
        return _detailedProcessingRecord;
    }

    public float[] getLatCorners() {
        if(_mapProjRecord == null) return null;

        final double latUL = _mapProjRecord.getAttributeDouble("1st line 1st pixel geodetic latitude");
        final double latUR = _mapProjRecord.getAttributeDouble("1st line last valid pixel geodetic latitude");
        final double latLL = _mapProjRecord.getAttributeDouble("Last line 1st pixel geodetic latitude");
        final double latLR = _mapProjRecord.getAttributeDouble("Last line last valid pixel geodetic latitude");
        return new float[]{(float)latUL, (float)latUR, (float)latLL, (float)latLR};
    }

    public float[] getLonCorners() {
        if(_mapProjRecord == null) return null;

        final double lonUL = _mapProjRecord.getAttributeDouble("1st line 1st pixel geodetic longitude");
        final double lonUR = _mapProjRecord.getAttributeDouble("1st line last valid pixel geodetic longitude");
        final double lonLL = _mapProjRecord.getAttributeDouble("Last line 1st pixel geodetic longitude");
        final double lonLR = _mapProjRecord.getAttributeDouble("Last line last valid pixel geodetic longitude");
        return new float[]{(float)lonUL, (float)lonUR, (float)lonLL, (float)lonLR};
    }

    public void addMetadata(MetadataElement sphElem) {

        CeosHelper.addMetadata(sphElem, _leaderFDR, "File Descriptor");
        CeosHelper.addMetadata(sphElem, _sceneHeaderRecord, "Scene Parameters");
        CeosHelper.addMetadata(sphElem, _mapProjRecord, "Map Projection");
        CeosHelper.addMetadata(sphElem, _platformPositionRecord, "Platform Position");
        CeosHelper.addMetadata(sphElem, _dataQualityRecord, "Data Quality");
        CeosHelper.addMetadata(sphElem, _histogramRecord, "Histogram");
        CeosHelper.addMetadata(sphElem, _attitudeRecord, "Attitude");
        CeosHelper.addMetadata(sphElem, _radiometricRecord, "Radiometric");
        CeosHelper.addMetadata(sphElem, _radiometricCompRecord, "Radiometric Compensation");
        CeosHelper.addMetadata(sphElem, _detailedProcessingRecord, "Detailed Processing");
        CeosHelper.addMetadata(sphElem, _facilityRecord, "Facility Related");
    }
}