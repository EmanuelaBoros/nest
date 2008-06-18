/*
 * $Id: FilePointerRecord.java,v 1.3 2008-06-17 20:35:10 lveci Exp $
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
package org.esa.nest.dataio.ceos.records;

import org.esa.nest.dataio.ceos.CeosFileReader;
import org.esa.nest.dataio.ceos.IllegalCeosFormatException;
import org.esa.beam.framework.datamodel.MetadataElement;

import java.io.IOException;

public class FilePointerRecord extends BaseRecord {

    public static final int RECORD_LENGTH = 360;

    private static final String IMAGE_FILE_CLASS_CODE = "IMGY";

  /*  private final String _codeCharacter;
    private final int _filePointerNumber;
    private final String _fileID;
    private final String _fileClass;
    private final String _fileClassCode;
    private final String _fileDataType;
    private final String _fileDataTypeCode;
    private final int _numberOfRecords;
    private final int _firstRecordLength;
    private final int _maxRecordLength;
    private final String _recordLengthType;
    private final String _recordLengthTypeCode;
    private final int _firstRecordVolumeNumber;
    private final int _finalRecordVolumeNumber;
    private final int _firstRecordNumberOfReferencedFile;       */

    private static String format = "ers";
    private static String recordDefinitionFile = "file_pointer_record.xml";

    public FilePointerRecord(final CeosFileReader reader) throws IOException, IllegalCeosFormatException {
        this(reader, -1);
    }

    public FilePointerRecord(final CeosFileReader reader, final long startPos) throws IOException,
                                                                                      IllegalCeosFormatException {
        super(reader, startPos, format, recordDefinitionFile);
        /*
        _codeCharacter = reader.readAn(2);
        reader.skipBytes(2);    // blank
        _filePointerNumber = reader.readI4();
        _fileID = reader.readAn(16);
        _fileClass = reader.readAn(28);
        _fileClassCode = reader.readAn(4);
        _fileDataType = reader.readAn(28);
        _fileDataTypeCode = reader.readAn(4);
        _numberOfRecords = (int) reader.readIn(8);
        _firstRecordLength = (int) reader.readIn(8);
        _maxRecordLength = (int) reader.readIn(8);
        _recordLengthType = reader.readAn(12);
        _recordLengthTypeCode = reader.readAn(4);
        _firstRecordVolumeNumber = (int) reader.readIn(2);
        _finalRecordVolumeNumber = (int) reader.readIn(2);
        _firstRecordNumberOfReferencedFile = (int) reader.readIn(8);
        // skip the last 208 blanks
        reader.skipBytes(208);  */
    }

  /*  public String getCodeCharacter() {
        return _codeCharacter;
    }

    public String getFileClass() {
        return _fileClass;
    }

    public String getFileClassCode() {
        return _fileClassCode;
    }

    public String getFileDataType() {
        return _fileDataType;
    }

    public String getFileDataTypeCode() {
        return _fileDataTypeCode;
    }

    public String getFileID() {
        return _fileID;
    }

    public int getFilePointerNumber() {
        return _filePointerNumber;
    }

    public int getFirstRecordLength() {
        return _firstRecordLength;
    }

    public int getMaxRecordLength() {
        return _maxRecordLength;
    }

    public int getNumberOfRecords() {
        return _numberOfRecords;
    }

    public String getRecordLengthType() {
        return _recordLengthType;
    }

    public String getRecordLengthTypeCode() {
        return _recordLengthTypeCode;
    }

    public int getFirstRecordVolumeNumber() {
        return _firstRecordVolumeNumber;
    }

    public int getFinalRecordVolumeNumber() {
        return _finalRecordVolumeNumber;
    }

    public int getFirstRecordNumberOfReferencedFile() {
        return _firstRecordNumberOfReferencedFile;
    }         */

    public boolean isImageFileRecord() {
        return FilePointerRecord.IMAGE_FILE_CLASS_CODE.equalsIgnoreCase(getAttributeString("File class code"));
    }

    public void assignMetadataTo(final MetadataElement root, final String suffix) {
        final MetadataElement elem = createMetadataElement("FilePointerRecord", suffix);
        root.addElement(elem);

        super.assignMetadataTo(elem);

    /*    elem.setAttributeString("Code character", _codeCharacter);
        elem.setAttributeInt("File pointer number", _filePointerNumber);
        elem.setAttributeString("File ID", _fileID);
        elem.setAttributeString("File class", _fileClass);
        elem.setAttributeString("File class code", _fileClassCode);
        elem.setAttributeString("File datatype", _fileDataType);
        elem.setAttributeString("File datatype code", _fileDataTypeCode);
        elem.setAttributeInt("Number of records", _numberOfRecords);
        elem.setAttributeInt("First record length", _firstRecordLength);
        elem.setAttributeInt("Max record length", _maxRecordLength);
        elem.setAttributeString("Record lengthtype", _recordLengthType);
        elem.setAttributeString("Record lengthtype code", _recordLengthTypeCode);
        elem.setAttributeInt("First record volume numer", _firstRecordVolumeNumber);
        elem.setAttributeInt("Final record volume number", _finalRecordVolumeNumber);
        elem.setAttributeInt("First record number of referenced file", _firstRecordNumberOfReferencedFile); */
    }
}
