/*
 * $Id: ProductDataUShortTest.java,v 1.3 2010-03-31 13:59:56 lveci Exp $
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
package org.esa.beam.framework.datamodel;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;
import org.esa.beam.GlobalTestConfig;
import org.esa.beam.util.SystemUtils;

import javax.imageio.stream.FileImageInputStream;
import javax.imageio.stream.FileImageOutputStream;
import java.io.File;
import java.io.IOException;

public class ProductDataUShortTest extends TestCase {

    private FileImageInputStream _inputStream;
    private FileImageOutputStream _outputStream;

    @Override
    protected void setUp() throws IOException {
        File outputFile = GlobalTestConfig.getBeamTestDataOutputFile("ProductData");
        outputFile.mkdirs();
        File streamFile = new File(outputFile, "ushort.img");
        streamFile.createNewFile();
        _inputStream = new FileImageInputStream(streamFile);
        _outputStream = new FileImageOutputStream(streamFile);
        assertNotNull(_inputStream);
        assertNotNull(_outputStream);
    }

    @Override
    protected void tearDown() {
        try {
            _inputStream.close();
            _outputStream.close();
        } catch (IOException e) {
        }
        SystemUtils.deleteFileTree(GlobalTestConfig.getBeamTestDataOutputDirectory());
    }

    public void testSingleValueConstructor() {
        ProductData instance = ProductData.createInstance(ProductData.TYPE_UINT16);
        instance.setElems(new short[]{-1});

        assertEquals(ProductData.TYPE_UINT16, instance.getType());
        assertEquals(65535, instance.getElemInt());
        assertEquals(65535L, instance.getElemUInt());
        assertEquals(65535.0F, instance.getElemFloat(), 0.0e-12F);
        assertEquals(65535.0D, instance.getElemDouble(), 0.0e-12D);
        assertEquals("65535", instance.getElemString());
        assertEquals(1, instance.getNumElems());
        Object data = instance.getElems();
        assertEquals(true, data instanceof short[]);
        assertEquals(1, ((short[]) data).length);
        assertEquals(true, instance.isScalar());
        assertEquals(true, instance.isInt());
        assertEquals("65535", instance.toString());

        ProductData expectedEqual = ProductData.createInstance(ProductData.TYPE_UINT16);
        expectedEqual.setElems(new short[]{-1});
        assertEquals(true, instance.equalElems(expectedEqual));

        ProductData expectedUnequal = ProductData.createInstance(ProductData.TYPE_UINT16);
        expectedUnequal.setElems(new short[]{-2});
        assertEquals(false, instance.equalElems(expectedUnequal));

//        StreamTest
        ProductData dataFromStream = null;
        try {
            instance.writeTo(_outputStream);
            dataFromStream = ProductData.createInstance(ProductData.TYPE_UINT16);
            dataFromStream.readFrom(_inputStream);
        } catch (IOException e) {
            fail("IOException not expected");
        }
        assertEquals(true, instance.equalElems(dataFromStream));
    }

    public void testConstructor() {
        ProductData instance = ProductData.createInstance(ProductData.TYPE_UINT16, 3);
        instance.setElems(new short[]{-1, 32767, -32768});

        assertEquals(ProductData.TYPE_UINT16, instance.getType());
        assertEquals(65535, instance.getElemIntAt(0));
        assertEquals(32767, instance.getElemIntAt(1));
        assertEquals(32768, instance.getElemIntAt(2));
        assertEquals(65535L, instance.getElemUIntAt(0));
        assertEquals(32767L, instance.getElemUIntAt(1));
        assertEquals(32768L, instance.getElemUIntAt(2));
        assertEquals(65535.0F, instance.getElemFloatAt(0), 0.0e-12F);
        assertEquals(32767.0F, instance.getElemFloatAt(1), 0.0e-12F);
        assertEquals(32768.0F, instance.getElemFloatAt(2), 0.0e-12F);
        assertEquals(65535.0D, instance.getElemDoubleAt(0), 0.0e-12D);
        assertEquals(32767.0D, instance.getElemDoubleAt(1), 0.0e-12D);
        assertEquals(32768.0D, instance.getElemDoubleAt(2), 0.0e-12D);
        assertEquals("65535", instance.getElemStringAt(0));
        assertEquals("32767", instance.getElemStringAt(1));
        assertEquals("32768", instance.getElemStringAt(2));
        assertEquals(3, instance.getNumElems());
        Object data2 = instance.getElems();
        assertEquals(true, data2 instanceof short[]);
        assertEquals(3, ((short[]) data2).length);
        assertEquals(false, instance.isScalar());
        assertEquals(true, instance.isInt());
        assertEquals("65535,32767,32768", instance.toString());

        ProductData expectedEqual = ProductData.createInstance(ProductData.TYPE_UINT16, 3);
        expectedEqual.setElems(new short[]{-1, 32767, -32768});
        assertEquals(true, instance.equalElems(expectedEqual));

        ProductData expectedUnequal = ProductData.createInstance(ProductData.TYPE_UINT16, 3);
        expectedUnequal.setElems(new short[]{-1, 32767, -32767});
        assertEquals(false, instance.equalElems(expectedUnequal));

//        StreamTest
        ProductData dataFromStream = null;
        try {
            instance.writeTo(_outputStream);
            dataFromStream = ProductData.createInstance(ProductData.TYPE_UINT16, 3);
            dataFromStream.readFrom(_inputStream);
        } catch (IOException e) {
            fail("IOException not expected");
        }
        assertEquals(true, instance.equalElems(dataFromStream));
    }

    public void testSetElemsAsString() {
        final ProductData pd = ProductData.createInstance(ProductData.TYPE_UINT16, 3);
        pd.setElems(new String[]{
                String.valueOf(Short.MAX_VALUE * 2 + 1),
                String.valueOf(Short.MAX_VALUE),
                String.valueOf(0),
        });

        assertEquals(Short.MAX_VALUE * 2 + 1, pd.getElemIntAt(0));
        assertEquals(Short.MAX_VALUE, pd.getElemIntAt(1));
        assertEquals(0, pd.getElemIntAt(2));
    }

    public void testSetElemsAsString_OutOfRange() {
        final ProductData pd1 = ProductData.createInstance(ProductData.TYPE_UINT16, 1);
        try {
            pd1.setElems(new String[]{String.valueOf(Short.MAX_VALUE * 2 + 2)});
        } catch (Exception e) {
            assertEquals(NumberFormatException.class, e.getClass());
            assertEquals("Value out of range. The value:'65536' is not an unsigned short value.", e.getMessage());
        }

        final ProductData pd2 = ProductData.createInstance(ProductData.TYPE_UINT16, 1);
        try {
            pd2.setElems(new String[]{String.valueOf(-1)});
        } catch (Exception e) {
            assertEquals(NumberFormatException.class, e.getClass());
            assertEquals("Value out of range. The value:'-1' is not an unsigned short value.", e.getMessage());
        }
    }
}