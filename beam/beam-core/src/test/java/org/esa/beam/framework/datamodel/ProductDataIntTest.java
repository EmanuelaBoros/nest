/*
 * $Id: ProductDataIntTest.java,v 1.1 2009-04-28 14:39:33 lveci Exp $
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

import java.io.File;
import java.io.IOException;

import javax.imageio.stream.FileImageInputStream;
import javax.imageio.stream.FileImageOutputStream;

import junit.framework.TestCase;

import org.esa.beam.GlobalTestConfig;
import org.esa.beam.util.SystemUtils;

public class ProductDataIntTest extends TestCase {

    private FileImageInputStream _inputStream;
    private FileImageOutputStream _outputStream;

    @Override
	protected void setUp() {
        File outputFile = GlobalTestConfig.getBeamTestDataOutputFile("ProductData");
        outputFile.mkdirs();
        File streamFile = new File(outputFile, "stream.img");
        try {
            streamFile.createNewFile();
            _inputStream = new FileImageInputStream(streamFile);
            _outputStream = new FileImageOutputStream(streamFile);
        } catch (IOException e) {
        }
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
        ProductData instance = ProductData.createInstance(ProductData.TYPE_INT32);
        instance.setElems(new int[]{2147483647});

        assertEquals(ProductData.TYPE_INT32, instance.getType());
        assertEquals(2147483647, instance.getElemInt());
        assertEquals(2147483647L, instance.getElemUInt());
        assertEquals(2147483647.0F, instance.getElemFloat(), 0.0e-12F);
        assertEquals(2147483647.0D, instance.getElemDouble(), 0.0e-12D);
        assertEquals("2147483647", instance.getElemString());
        assertEquals(true, instance.getElemBoolean());
        assertEquals(1, instance.getNumElems());
        Object data = instance.getElems();
        assertEquals(true, data instanceof int[]);
        assertEquals(1, ((int[]) data).length);
        assertEquals(true, instance.isScalar());
        assertEquals(true, instance.isInt());
        assertEquals("2147483647", instance.toString());

        ProductData expectedEqual = ProductData.createInstance(ProductData.TYPE_INT32);
        expectedEqual.setElems(new int[]{2147483647});
        assertEquals(true, instance.equalElems(expectedEqual));

        ProductData expectedUnequal = ProductData.createInstance(ProductData.TYPE_INT32);
        expectedUnequal.setElems(new int[]{2147483646});
        assertEquals(false, instance.equalElems(expectedUnequal));

//        StreamTest
        ProductData dataFromStream = null;
        try {
            instance.writeTo(_outputStream);
            dataFromStream = ProductData.createInstance(ProductData.TYPE_INT32);
            dataFromStream.readFrom(_inputStream);
        } catch (IOException e) {
            fail("IOException not expected");
        }
        assertEquals(true, instance.equalElems(dataFromStream));
    }

    public void testConstructor() {
        ProductData instance = ProductData.createInstance(ProductData.TYPE_INT32, 3);
        instance.setElems(new int[]{-1, 2147483647, -2147483648});

        assertEquals(ProductData.TYPE_INT32, instance.getType());
        assertEquals(-1, instance.getElemIntAt(0));
        assertEquals(2147483647, instance.getElemIntAt(1));
        assertEquals(-2147483648, instance.getElemIntAt(2));
        assertEquals(-1L, instance.getElemUIntAt(0));
        assertEquals(2147483647L, instance.getElemUIntAt(1));
        assertEquals(-2147483648L, instance.getElemUIntAt(2));
        assertEquals(-1.0F, instance.getElemFloatAt(0), 0.0e-12F);
        assertEquals(2147483647.0F, instance.getElemFloatAt(1), 0.0e-12F);
        assertEquals(-2147483648.0F, instance.getElemFloatAt(2), 0.0e-12F);
        assertEquals(-1.0D, instance.getElemDoubleAt(0), 0.0e-12D);
        assertEquals(2147483647.0D, instance.getElemDoubleAt(1), 0.0e-12D);
        assertEquals(-2147483648.0D, instance.getElemDoubleAt(2), 0.0e-12D);
        assertEquals("-1", instance.getElemStringAt(0));
        assertEquals("2147483647", instance.getElemStringAt(1));
        assertEquals("-2147483648", instance.getElemStringAt(2));
        assertEquals(true, instance.getElemBooleanAt(0));
        assertEquals(true, instance.getElemBooleanAt(1));
        assertEquals(true, instance.getElemBooleanAt(2));
        assertEquals(3, instance.getNumElems());
        Object data2 = instance.getElems();
        assertEquals(true, data2 instanceof int[]);
        assertEquals(3, ((int[]) data2).length);
        assertEquals(false, instance.isScalar());
        assertEquals(true, instance.isInt());
        assertEquals("-1,2147483647,-2147483648", instance.toString());

        ProductData expectedEqual = ProductData.createInstance(ProductData.TYPE_INT32, 3);
        expectedEqual.setElems(new int[]{-1, 2147483647, -2147483648});
        assertEquals(true, instance.equalElems(expectedEqual));

        ProductData expectedUnequal = ProductData.createInstance(ProductData.TYPE_INT32, 3);
        expectedUnequal.setElems(new int[]{-1, 2147483647, -2147483647});
        assertEquals(false, instance.equalElems(expectedUnequal));

//        StreamTest
        ProductData dataFromStream = null;
        try {
            instance.writeTo(_outputStream);
            dataFromStream = ProductData.createInstance(ProductData.TYPE_INT32, 3);
            dataFromStream.readFrom(_inputStream);
        } catch (IOException e) {
            fail("IOException not expected");
        }
        assertEquals(true, instance.equalElems(dataFromStream));
    }

    public void testSetElemsAsString() {
        final ProductData pd = ProductData.createInstance(ProductData.TYPE_INT32, 3);
        pd.setElems(new String[]{
                String.valueOf(Integer.MAX_VALUE),
                String.valueOf(Integer.MIN_VALUE),
                String.valueOf(0),
        });
        
        assertEquals(Integer.MAX_VALUE, pd.getElemIntAt(0));
        assertEquals(Integer.MIN_VALUE, pd.getElemIntAt(1));
        assertEquals(0, pd.getElemIntAt(2));
}
    public void testSetElemsAsString_OutOfRange() {
        final ProductData pd1 = ProductData.createInstance(ProductData.TYPE_INT32, 1);
        try {
            pd1.setElems(new String[]{String.valueOf((long)Integer.MAX_VALUE + 1)});
        } catch (Exception e) {
            assertEquals(NumberFormatException.class, e.getClass());
            assertEquals(true, e.getMessage().contains("2147483648"));
        }

        final ProductData pd2 = ProductData.createInstance(ProductData.TYPE_INT32, 1);
        try {
            pd2.setElems(new String[]{String.valueOf((long)Integer.MIN_VALUE -1)});
        } catch (Exception e) {
            assertEquals(NumberFormatException.class, e.getClass());
            assertEquals(true, e.getMessage().contains("-2147483649"));
        }
    }
}
