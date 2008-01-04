/*
 * $Id: Ancillary2Record.java,v 1.1 2008-01-04 16:23:10 lveci Exp $
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

import java.io.IOException;

public class Ancillary2Record extends BaseRecord {

    private String _sensorOperationMode;
    private int _lowerLimitOfStrengthAfterCorrection;
    private int _upperLimitOfStrengthAfterCorrection;
    private String _sensorGains;

    public Ancillary2Record(final CeosFileReader reader) throws IOException, IllegalCeosFormatException {
        this(reader, -1);
    }

    public Ancillary2Record(final CeosFileReader reader, final long startPos) throws IOException,
                                                                                     IllegalCeosFormatException {
        super(reader, startPos);

        readGeneralFields(reader);

        reader.seek(getAbsolutPosition(getRecordLength()));
    }

    private void readGeneralFields(final CeosFileReader reader) throws IOException,
                                                                       IllegalCeosFormatException {
        _sensorOperationMode = reader.readAn(4);
        _lowerLimitOfStrengthAfterCorrection = reader.readI4();
        _upperLimitOfStrengthAfterCorrection = reader.readI4();
        reader.skipBytes(32);   // skip 30 + 1 + 1
        _sensorGains = reader.readAn(6);

        readSpecificFields(reader);
    }

    protected void readSpecificFields(final CeosFileReader reader) throws IOException,
                                                                          IllegalCeosFormatException {
    }

    public String getSensorOperationMode() {
        return _sensorOperationMode;
    }

    public int getLowerLimitOfStrengthAfterCorrection() {
        return _lowerLimitOfStrengthAfterCorrection;
    }

    public int getUpperLimitOfStrengthAfterCorrection() {
        return _upperLimitOfStrengthAfterCorrection;
    }

    public String getSensorGains() {
        return _sensorGains;
    }
}
