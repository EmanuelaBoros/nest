/*
 * $Id: Ancillary3Record.java,v 1.2 2008-01-07 15:04:28 lveci Exp $
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

public class Ancillary3Record extends BaseRecord {

    private final String _orbitalElementsType;
    private final int _numDataPoints;
    private final int _firstPointYear;
    private final int _firstPointMonth;
    private final int _firstPointDay;
    private final int _firstPointTotalDays;
    private final double _firstPointTotalSeconds;
    private final double _intervalTimeBetweenPoints;
    private final String _referenceCoordinateSystem;
    private final double _positionalErrorFlightDirection;
    private final double _positionalErrorFlightVerticalDirection;
    private final double _positionalErrorRadiusDirection;
    private final double _velocityErrorFlightDirection;
    private final double _velocityErrorFlightVerticalDirection;
    private final double _velocityErrorRadiusDirection;
    private final DataPoint[] _dataPoints;
    private final int _flagLeapSecond;

    public Ancillary3Record(final CeosFileReader reader) throws IOException,
                                                                IllegalCeosFormatException {
        this(reader, -1);
    }

    public Ancillary3Record(final CeosFileReader reader, final long startPos) throws IOException,
                                                                                     IllegalCeosFormatException {
        super(reader, startPos);

        _orbitalElementsType = reader.readAn(32);
        reader.skipBytes(16 * 6); // 6 x orbital element [BLANK]
        _numDataPoints = reader.readI4();
        _firstPointYear = reader.readI4();

        //tmp
            _firstPointMonth = 0;
            _firstPointDay = 0;
            _firstPointTotalDays = 0;
            _firstPointTotalSeconds = 0;
            _intervalTimeBetweenPoints =0;
            _referenceCoordinateSystem = "";
            _positionalErrorFlightDirection = 0;
            _positionalErrorFlightVerticalDirection = 0;
            _positionalErrorRadiusDirection = 0;
            _velocityErrorFlightDirection =0;
            _velocityErrorFlightVerticalDirection = 0;
            _velocityErrorRadiusDirection =0;
           _flagLeapSecond = 0;
           _dataPoints = new DataPoint[28];

     /*   _firstPointMonth = reader.readI4();
        _firstPointDay = reader.readI4();
        _firstPointTotalDays = reader.readI4();
        _firstPointTotalSeconds = reader.readFn(22);
        _intervalTimeBetweenPoints = reader.readFn(22);
        _referenceCoordinateSystem = reader.readAn(64);
        reader.skipBytes(22); // greenwichMeanHourAngle [BLANK]
        _positionalErrorFlightDirection = reader.readFn(16);
        _positionalErrorFlightVerticalDirection = reader.readFn(16);
        _positionalErrorRadiusDirection = reader.readFn(16);
        _velocityErrorFlightDirection = reader.readFn(16);
        _velocityErrorFlightVerticalDirection = reader.readFn(16);
        _velocityErrorRadiusDirection = reader.readFn(16);
        _dataPoints = new DataPoint[28];
        for (int i = 0; i < _dataPoints.length; i++) {
            _dataPoints[i] = new DataPoint(reader.readFn(22), reader.readFn(22), reader.readFn(22),
                                           reader.readFn(22), reader.readFn(22), reader.readFn(22));
        }
        reader.skipBytes(18);
        _flagLeapSecond = (int) reader.readIn(1);
        reader.skipBytes(579);      */
    }

    public String getOrbitalElementsType() {
        return _orbitalElementsType;
    }

//    public double getOrbitalElement1() {
//        return _orbitalElement1;
//    }
//
//    public double getOrbitalElement2() {
//        return _orbitalElement2;
//    }
//
//    public double getOrbitalElement3() {
//        return _orbitalElement3;
//    }
//
//    public double getOrbitalElement4() {
//        return _orbitalElement4;
//    }
//
//    public double getOrbitalElement5() {
//        return _orbitalElement5;
//    }
//
//    public double getOrbitalElement6() {
//        return _orbitalElement6;
//    }

    public int getNumDataPoints() {
        return _numDataPoints;
    }

    public int getFirstPointYear() {
        return _firstPointYear;
    }

    public int getFirstPointMonth() {
        return _firstPointMonth;
    }

    public int getFirstPointDay() {
        return _firstPointDay;
    }

    public int getFirstPointTotalDays() {
        return _firstPointTotalDays;
    }

    public double getFirstPointTotalSeconds() {
        return _firstPointTotalSeconds;
    }

    public double getIntervalTimeBetweenPoints() {
        return _intervalTimeBetweenPoints;
    }

    public String getReferenceCoordinateSystem() {
        return _referenceCoordinateSystem;
    }

    public double getPositionalErrorFlightDirection() {
        return _positionalErrorFlightDirection;
    }

    public double getPositionalErrorFlightVerticalDirection() {
        return _positionalErrorFlightVerticalDirection;
    }

    public double getPositionalErrorRadiusDirection() {
        return _positionalErrorRadiusDirection;
    }

    public double getVelocityErrorFlightDirection() {
        return _velocityErrorFlightDirection;
    }

    public double getVelocityErrorFlightVerticalDirection() {
        return _velocityErrorFlightVerticalDirection;
    }

    public double getVelocityErrorRadiusDirection() {
        return _velocityErrorRadiusDirection;
    }

    public DataPoint[] getDataPoints() {
        return _dataPoints;
    }

    public int getFlagLeapSecond() {
        return _flagLeapSecond;
    }

    public static class DataPoint {

        private final double _positionalVectorDataPointX;
        private final double _positionalVectorDataPointY;
        private final double _positionalVectorDataPointZ;
        private final double _velocityVectorDataPointX;
        private final double _velocityVectorDataPointY;
        private final double _velocityVectorDataPointZ;

        public DataPoint(final double positionalVectorDataPointX, final double positionalVectorDataPointY,
                         final double positionalVectorDataPointZ,
                         final double velocityVectorDataPointX, final double velocityVectorDataPointY,
                         final double velocityVectorDataPointZ) {
            _positionalVectorDataPointX = positionalVectorDataPointX;
            _positionalVectorDataPointY = positionalVectorDataPointY;
            _positionalVectorDataPointZ = positionalVectorDataPointZ;
            _velocityVectorDataPointX = velocityVectorDataPointX;
            _velocityVectorDataPointY = velocityVectorDataPointY;
            _velocityVectorDataPointZ = velocityVectorDataPointZ;
        }

        public double getPositionalVectorDataPointX() {
            return _positionalVectorDataPointX;
        }

        public double getPositionalVectorDataPointY() {
            return _positionalVectorDataPointY;
        }

        public double getPositionalVectorDataPointZ() {
            return _positionalVectorDataPointZ;
        }

        public double getVelocityVectorDataPointX() {
            return _velocityVectorDataPointX;
        }

        public double getVelocityVectorDataPointY() {
            return _velocityVectorDataPointY;
        }

        public double getVelocityVectorDataPointZ() {
            return _velocityVectorDataPointZ;
        }
    }
}
