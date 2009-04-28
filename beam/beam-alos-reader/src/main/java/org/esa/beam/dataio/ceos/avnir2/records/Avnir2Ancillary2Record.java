package org.esa.beam.dataio.ceos.avnir2.records;

import org.esa.beam.dataio.ceos.CeosFileReader;
import org.esa.beam.dataio.ceos.IllegalCeosFormatException;
import org.esa.beam.dataio.ceos.records.Ancillary2Record;

import java.io.IOException;

/*
 * $Id: Avnir2Ancillary2Record.java,v 1.1 2009-04-28 14:37:13 lveci Exp $
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
 * Created by Marco Peters.
 *
 * @author Marco Peters
 * @version $Revision: 1.1 $ $Date: 2009-04-28 14:37:13 $
 */
public class Avnir2Ancillary2Record extends Ancillary2Record {

    private long[] _bandExposureCoefficients;
    private double[] _bandDetectorTemperature;
    private double[] _bandDetectorAssemblyTemperature;
    private double[] _bandGains;
    private double[] _bandOffsets;
    private double _signalProcessingUnitTemperature;

    public Avnir2Ancillary2Record(final CeosFileReader reader) throws IOException,
                                                                      IllegalCeosFormatException {
        this(reader, -1);
    }

    public Avnir2Ancillary2Record(final CeosFileReader reader, final long startPos) throws IOException,
                                                                                           IllegalCeosFormatException {
        super(reader, startPos);
    }

    @Override
    protected void readSpecificFields(final CeosFileReader reader) throws IOException,
                                                                          IllegalCeosFormatException {
        reader.seek(getAbsolutPosition(24));
        _bandExposureCoefficients = new long[4];
        _bandExposureCoefficients[0] = reader.readIn(5);
        _bandExposureCoefficients[1] = reader.readIn(5);
        _bandExposureCoefficients[2] = reader.readIn(5);
        _bandExposureCoefficients[3] = reader.readIn(5);

        reader.seek(getAbsolutPosition(78));
        _bandDetectorTemperature = new double[4];
        _bandDetectorTemperature[0] = reader.readFn(8);
        _bandDetectorTemperature[1] = reader.readFn(8);
        _bandDetectorTemperature[2] = reader.readFn(8);
        _bandDetectorTemperature[3] = reader.readFn(8);

        _bandDetectorAssemblyTemperature = new double[4];
        _bandDetectorAssemblyTemperature[0] = reader.readFn(8);
        _bandDetectorAssemblyTemperature[1] = reader.readFn(8);
        _bandDetectorAssemblyTemperature[2] = reader.readFn(8);
        _bandDetectorAssemblyTemperature[3] = reader.readFn(8);

        _signalProcessingUnitTemperature = reader.readFn(8);

        reader.seek(getAbsolutPosition(2702));
        _bandGains = new double[4];
        _bandOffsets = new double[4];
        _bandGains[0] = reader.readFn(8);
        _bandOffsets[0] = reader.readFn(8);
        _bandGains[1] = reader.readFn(8);
        _bandOffsets[1] = reader.readFn(8);
        _bandGains[2] = reader.readFn(8);
        _bandOffsets[2] = reader.readFn(8);
        _bandGains[3] = reader.readFn(8);
        _bandOffsets[3] = reader.readFn(8);
    }

    public double getDetectorAssemblyTemperature(final int bandNumber) {
        return _bandDetectorAssemblyTemperature[bandNumber - 1];
    }

    public double getDetectorTemperature(final int bandNumber) {
        return _bandDetectorTemperature[bandNumber - 1];
    }

    public long getExposureCoefficient(final int bandNumber) {
        return _bandExposureCoefficients[bandNumber - 1];
    }

    public double getSignalProcessingUnitTemperature() {
        return _signalProcessingUnitTemperature;
    }

    public double getAbsoluteCalibrationGain(final int bandNumber) {
        return _bandGains[bandNumber - 1];
    }

    public double getAbsoluteCalibrationOffset(final int bandNumber) {
        return _bandOffsets[bandNumber - 1];
    }
}
