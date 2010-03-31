package org.esa.beam.dataio.ceos.records;

/*
 * $Id: BaseTrailerRecordTest.java,v 1.2 2010-03-31 13:59:55 lveci Exp $
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

import junit.framework.TestCase;
import org.esa.beam.dataio.ceos.CeosFileReader;
import org.esa.beam.dataio.ceos.CeosTestHelper;
import org.esa.beam.dataio.ceos.IllegalCeosFormatException;

import javax.imageio.stream.ImageOutputStream;
import javax.imageio.stream.MemoryCacheImageOutputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;

/**
 * Created by Marco.
 *
 * @author Marco
 * @version $Revision: 1.2 $ $Date: 2010-03-31 13:59:55 $
 */
public abstract class BaseTrailerRecordTest extends TestCase {

    private String _prefix;
    private CeosFileReader _reader;

    @Override
    protected void setUp() throws Exception {
        final ByteArrayOutputStream os = new ByteArrayOutputStream(24);
        MemoryCacheImageOutputStream ios = new MemoryCacheImageOutputStream(os);
        _prefix = "TrailerRecordTest_prefix";
        ios.writeBytes(_prefix);
        writeRecordData(ios);
        ios.writeBytes("TrailerRecordTest_suffix"); // as suffix
        _reader = new CeosFileReader(ios);
    }

    public void testInit_SimpleConstructor() throws IOException,
                                                    IllegalCeosFormatException {
        _reader.seek(_prefix.length());

        final BaseTrailerRecord record = createTrailerRecord(_reader);

        assertRecord(record);
    }


    public void testInit() throws IOException,
                                  IllegalCeosFormatException {
        final BaseTrailerRecord record = createTrailerRecord(_reader, _prefix.length());

        assertRecord(record);
    }

    private void writeRecordData(final ImageOutputStream ios) throws IOException {
        BaseRecordTest.writeRecordData(ios);

        ios.writeBytes("   1"); // number of trailer records // I4
        ios.writeBytes("   1"); // number of trailer in one CCD // I4

        writeHistograms(ios);

        CeosTestHelper.writeBlanks(ios, 248);
    }

    private void assertRecord(final BaseTrailerRecord record) throws IOException {
        BaseRecordTest.assertRecord(record);
        assertEquals(_prefix.length(), record.getStartPos());
//        assertEquals(_prefix.length() + 8460, _ios.getStreamPosition());

        assertEquals(1, record.getNumTrailerRecords());
        assertEquals(1, record.getNumTrailerRecordsInOneCCDUnit());
        assertHistograms(record);
    }

    protected abstract BaseTrailerRecord createTrailerRecord(CeosFileReader reader) throws IOException,
                                                                                           IllegalCeosFormatException;

    protected abstract BaseTrailerRecord createTrailerRecord(final CeosFileReader reader, final int startPos) throws
                                                                                                              IOException,
                                                                                                              IllegalCeosFormatException;

    protected abstract void writeHistograms(final ImageOutputStream ios) throws IOException;

    protected abstract void assertHistograms(final BaseTrailerRecord record);

}