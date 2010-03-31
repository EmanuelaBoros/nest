/*
 * $Id: ModisUint8BandReader.java,v 1.2 2010-03-31 13:59:56 lveci Exp $
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
package org.esa.beam.dataio.modis.bandreader;

import ncsa.hdf.hdflib.HDFException;
import org.esa.beam.dataio.modis.hdf.lib.HDF;
import org.esa.beam.framework.datamodel.ProductData;

public class ModisUint8BandReader extends ModisBandReader {

    private byte[] _line;
    private short min;
    private short max;
    private byte fill;
    private byte[] targetData;
    private int targetIdx;

    public ModisUint8BandReader(final int sdsId, final int layer, final boolean is3d) {
        super(sdsId, layer, is3d);
    }

    /**
     * Retrieves the data type of the band
     *
     * @return always {@link org.esa.beam.framework.datamodel.ProductData#TYPE_UINT16}
     */
    @Override
    public int getDataType() {
        return ProductData.TYPE_UINT8;
    }

    @Override
    protected void prepareForReading(final int sourceOffsetX, final int sourceOffsetY, final int sourceWidth,
                                     final int sourceHeight, final int sourceStepX, final int sourceStepY,
                                     final ProductData destBuffer) {
        fill = (byte) Math.floor(_fillValue + 0.5);
        if (_validRange == null) {
            min = 0;
            max = Byte.MAX_VALUE * 2 + 1;
        } else {
            min = (short) Math.floor(_validRange.getMin() + 0.5);
            max = (short) Math.floor(_validRange.getMax() + 0.5);
        }
        targetData = (byte[]) destBuffer.getElems();
        targetIdx = 0;
        ensureLineWidth(sourceWidth);
    }

    @Override
    protected void readLine() throws HDFException {
        HDF.getWrap().SDreaddata(_sdsId, _start, _stride, _count, _line);
    }

    @Override
    protected void validate(final int x) {
        final int value = _line[x] & 0xff;
        if (value < min || value > max) {
            _line[x] = fill;
        }
    }

    @Override
    protected void assign(final int x) {
        targetData[targetIdx++] = _line[x];
    }

    private void ensureLineWidth(final int sourceWidth) {
        if ((_line == null) || (_line.length != sourceWidth)) {
            _line = new byte[sourceWidth];
        }
    }
}
