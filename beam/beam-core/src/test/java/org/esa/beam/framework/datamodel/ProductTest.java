/*
 * $Id: ProductTest.java,v 1.1 2009-04-28 14:39:33 lveci Exp $
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

import com.bc.ceres.core.ProgressMonitor;
import com.bc.jexp.ParseException;
import com.bc.jexp.Term;
import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;
import org.esa.beam.framework.dataio.AbstractProductReader;
import org.esa.beam.framework.dataio.DecodeQualification;
import org.esa.beam.framework.dataio.ProductReader;
import org.esa.beam.framework.dataio.ProductReaderPlugIn;
import org.esa.beam.framework.dataop.maptransf.*;
import org.esa.beam.util.BeamConstants;
import org.esa.beam.util.BitRaster;
import org.esa.beam.util.ObjectUtils;
import org.esa.beam.util.ProductUtilsTest;
import org.esa.beam.util.io.BeamFileFilter;

import java.awt.Dimension;
import java.io.IOException;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.Locale;

public class ProductTest extends TestCase {

    private static final String _prodType = "TestProduct";
    private static final int _sceneWidth = 20;
    private static final int _sceneHeight = 30;

    private Product _product;

//    public static void main(String[] args) {
//        junit.swingui.TestRunner.run(ProductTest.class);
//    }

    public ProductTest(String testName) {
        super(testName);
    }

    public static Test suite() {
        return new TestSuite(ProductTest.class);
    }

    /**
     * Initializes the tests
     */
    @Override
    protected void setUp() {
        _product = new Product("product", _prodType, _sceneWidth, _sceneHeight);
        _product.setModified(false);
    }

    public void testAcceptVisitor() {
        LinkedListProductVisitor visitor = new LinkedListProductVisitor();
        List<String> expectedList = new LinkedList<String>();
        assertEquals(expectedList, visitor.getVisitedList());

        _product.acceptVisitor(visitor);
        expectedList.add("flag_codings");
        expectedList.add("index_codings");
        expectedList.add("metadata");
        expectedList.add("product");
        assertEquals(expectedList, visitor.getVisitedList());

        try {
            _product.acceptVisitor(null);
            fail("Null argument for visitor not allowed");
        } catch (IllegalArgumentException e) {
        }
    }

    public void testAddProductNodeListener() {
        ProductNodeListener listener = new DummyProductNodeListener();
        assertEquals("addProductNodeListener null", false, _product.addProductNodeListener(null));
        assertEquals("addProductNodeListener", true, _product.addProductNodeListener(listener));
        assertEquals("addProductNodeListener contained listener", false, _product.addProductNodeListener(listener));
    }

    public void testFireNodeAdded() {
    }

    public void testFireNodeChanged() {
    }

    public void testFireNodeRemoved() {
    }

    public void testHasProductListeners() {
    }

    public void testRemoveProductListener() {
    }

    public void testSetWriter() {
    }

    public void testSetAndGetReader() {
        Product product = new Product("name", BeamConstants.MERIS_RR_L1B_PRODUCT_TYPE_NAME, 312, 213);

        assertNull(product.getProductReader());

        DummyProductReader reader1 = new DummyProductReader(new DummyProductReaderPlugIn());
        product.setProductReader(reader1);
        assertSame(reader1, product.getProductReader());

        DummyProductReader reader2 = new DummyProductReader(new DummyProductReaderPlugIn());
        product.setProductReader(reader2);
        assertSame(reader2, product.getProductReader());

        try {
            product.setProductReader(null);
            fail("IllegalArgumentException expected since the parameter is null");
        } catch (IllegalArgumentException e) {
            //IllegalArgumentException expected since the parameter is null
        }
    }

    public void testGetWriter() {
    }

    public void testAddBandWithBandParameters() {
        assertEquals(0, _product.getNumBands());
        assertEquals(0, _product.getBandNames().length);
        assertNull(_product.getBand("band1"));
        assertEquals(false, _product.containsBand("band1"));

        _product.addBand(new Band("band1", ProductData.TYPE_FLOAT32, _sceneWidth, _sceneHeight));

        assertEquals(1, _product.getNumBands());
        assertEquals("band1", _product.getBandNames()[0]);
        assertEquals(true, _product.containsBand("band1"));
        assertNotNull(_product.getBandAt(0));
        assertEquals("band1", _product.getBandAt(0).getName());
        assertNotNull(_product.getBand("band1"));
        assertEquals("band1", _product.getBand("band1").getName());
    }

    /**
     * Tests the functionality of getType
     */
    public void testGetType() {
        Product prod;

        prod = new Product("TestName", "TEST", _sceneWidth, _sceneHeight);
        assertEquals("TEST", prod.getProductType());

        prod = new Product("TestName", "TEST", _sceneWidth, _sceneHeight, null);
        assertEquals("TEST", prod.getProductType());
    }

    /**
     * Tests the functionality of getBandOutputRasterWidth
     */
    public void testGetSceneWidth() {
        Product prod;

        prod = new Product("TestName", _prodType, 243, _sceneHeight);
        assertEquals(243, prod.getSceneRasterWidth());

        prod = new Product("TestName", _prodType, 789, _sceneHeight, null);
        assertEquals(789, prod.getSceneRasterWidth());
    }

    /**
     * Tests the functionality of getBandOutputRasterHeight
     */
    public void testGetSceneHeight() {
        Product prod;

        prod = new Product("TestName", _prodType, _sceneWidth, 373);
        assertEquals(373, prod.getSceneRasterHeight());

        prod = new Product("TestName", _prodType, _sceneWidth, 427, null);
        assertEquals(427, prod.getSceneRasterHeight());
    }

    public void testBitmaskHandlingByte() {
        Product product = new Product("Y", "X", 4, 4);

        Band band = product.addBand("flags", ProductData.TYPE_INT8);
        band.setSynthetic(true);
        final byte F1 = 0x01;
        final byte F2 = 0x02;
        final byte F3 = 0x04;

        FlagCoding flagCoding = new FlagCoding("flags");
        flagCoding.addFlag("F1", F1, null);
        flagCoding.addFlag("F2", F2, null);
        flagCoding.addFlag("F3", F3, null);

        product.addFlagCoding(flagCoding);
        band.setFlagCoding(flagCoding);

        band.ensureRasterData();
        final byte[] elems = new byte[]{
                0, F1, F2, F3,
                F1, 0, F1 + F2, F1 + F3,
                F2, F1 + F2, 0, F2 + F3,
                F3, F1 + F3, F2 + F3, 0,
        };
        band.getRasterData().setElems(elems);
        product.setModified(false);

        final boolean[] F1_MASK = new boolean[]{
                false, true, false, false,
                true, false, true, true,
                false, true, false, false,
                false, true, false, false
        };
        testBitmaskHandling(product, "flags.F1", F1_MASK);

        final boolean[] F1_AND_F2_MASK = new boolean[]{
                false, false, false, false,
                false, false, true, false,
                false, true, false, false,
                false, false, false, false
        };
        testBitmaskHandling(product, "flags.F1 AND flags.F2", F1_AND_F2_MASK);

        final boolean[] F1_AND_F2_OR_F3_MASK = new boolean[]{
                false, false, false, true,
                false, false, true, true,
                false, true, false, true,
                true, true, true, false
        };
        testBitmaskHandling(product, "(flags.F1 AND flags.F2) OR flags.F3", F1_AND_F2_OR_F3_MASK);
    }

    public void testBitmaskHandlingUShort() {
        Product product = new Product("Y", "X", 4, 4);

        Band band = new Band("flags", ProductData.TYPE_UINT16, 4, 4);
        band.setSynthetic(true);
        product.addBand(band);

        final byte F1 = 0x0001;
        final byte F2 = 0x0002;
        final byte F3 = 0x0004;

        FlagCoding flagCoding = new FlagCoding("flags");
        flagCoding.addFlag("F1", F1, null);
        flagCoding.addFlag("F2", F2, null);
        flagCoding.addFlag("F3", F3, null);

        product.addFlagCoding(flagCoding);
        band.setFlagCoding(flagCoding);

        ProductData data = band.createCompatibleRasterData();
        final short[] elems = new short[]{
                0, F1, F2, F3,
                F1, 0, F1 + F2, F1 + F3,
                F2, F1 + F2, 0, F2 + F3,
                F3, F1 + F3, F2 + F3, 0
        };
        data.setElems(elems);
        band.setRasterData(data);
        product.setModified(false);

        final boolean[] F1_MASK = new boolean[]{
                false, true, false, false,
                true, false, true, true,
                false, true, false, false,
                false, true, false, false
        };
        testBitmaskHandling(product, "flags.F1", F1_MASK);

        final boolean[] F1_AND_F2_MASK = new boolean[]{
                false, false, false, false,
                false, false, true, false,
                false, true, false, false,
                false, false, false, false
        };
        testBitmaskHandling(product, "flags.F1 AND flags.F2", F1_AND_F2_MASK);

        final boolean[] F1_AND_F2_OR_F3_MASK = new boolean[]{
                false, false, false, true,
                false, false, true, true,
                false, true, false, true,
                true, true, true, false
        };
        testBitmaskHandling(product, "(flags.F1 AND flags.F2) OR flags.F3", F1_AND_F2_OR_F3_MASK);
    }

    private void testBitmaskHandling(Product product, String expr, boolean[] expected) {
        testBitmaskHandlingFullRaster(product, expr, expected);
        testBitmaskHandlingLineWise(product, expr, expected);
    }

    private static void testBitmaskHandlingFullRaster(Product product, String expr, boolean[] expected) {

        boolean[] res = new boolean[4 * 4];
        Term term = null;

        try {
            term = product.createTerm(expr);
        } catch (ParseException e) {
            fail("unexpected BitmaskExpressionParseException: " + e.getMessage());
        }

        try {
            product.readBitmask(0, 0, 4, 4, term, res, ProgressMonitor.NULL);
        } catch (IOException e) {
            fail("unexpected IOException: " + e.getMessage());
        }

        assertEquals(true, Arrays.equals(expected, res));
    }

    private void testBitmaskHandlingLineWise(Product product, String expr, boolean[] expected) {

        boolean[] res = new boolean[4 * 4];
        Term term = null;

        try {
            term = product.createTerm(expr);
        } catch (ParseException e) {
            fail("unexpected BitmaskExpressionParseException: " + e.getMessage());
        }

        try {
            readBitmaskLineWise(product, term, res);
        } catch (IOException e) {
            fail("unexpected IOException: " + e.getMessage());
        }

        assertEquals(true, Arrays.equals(expected, res));
    }

    private static void readBitmaskLineWise(Product product, Term term, boolean[] res) throws IOException {
        boolean[] line = new boolean[4];
        for (int y = 0; y < 4; y++) {
            product.readBitmask(0, y, 4, 1, term, line, ProgressMonitor.NULL);
            for (int x = 0; x < 4; x++) {
                res[y * 4 + x] = line[x];
            }
        }
    }

    public void testGetAndSetRefNo() {
        assertEquals(0, _product.getRefNo());

        try {
            _product.setRefNo(0);
            fail("IllegalArgumentException expected");
        } catch (IllegalArgumentException e) {
            // expectet if value out of range
        } catch (IllegalStateException e) {
            fail("IllegalStateException not expected");
        }

        try {
            _product.setRefNo(14);
        } catch (IllegalArgumentException e) {
            fail("IllegalArgumentException not expected");
        } catch (IllegalStateException e) {
            fail("IllegalStateException not expected");
        }
        assertEquals(14, _product.getRefNo());

        try {
            _product.setRefNo(23);
            fail("IllegalStateException expected");
        } catch (IllegalArgumentException e) {
            fail("IllegalArgumentException not expected");
        } catch (IllegalStateException e) {
            // expected if the reference number was alredy set
        }

        // no exception expected when the reference number to be set is the same as the one already set
        try {
            _product.setRefNo(14);
        } catch (IllegalArgumentException e) {
            fail("IllegalArgumentException not expected");
        } catch (IllegalStateException e) {
            fail("IllegalStateException not expected");
        }
    }


    public void testModifiedProperty() {

        assertEquals("product should be initially un-modified", false, _product.isModified());
        _product.setModified(true);
        assertEquals(true, _product.isModified());
        _product.setModified(false);
        assertEquals(false, _product.isModified());
    }

    public void testModifiedFlagAfterBandHasBeenAddedAndRemoved() {

        Band band = new Band("band1", ProductData.TYPE_FLOAT32, _sceneWidth, _sceneHeight);

        assertEquals(null, _product.getBand("band1"));

        //
        _product.addBand(band);
        assertEquals(band, _product.getBand("band1"));
        assertEquals("added band, modified flag should be set", true, _product.isModified());
        _product.setModified(false);

        _product.removeBand(band);
        assertEquals(null, _product.getBand("band1"));
        assertEquals("removed band, modified flag should be set", true, _product.isModified());
    }

    public void testModifiedFlagAfterBandHasBeenModified() {

        Band band = new Band("band1", ProductData.TYPE_FLOAT32, _sceneWidth, _sceneHeight);
        _product.addBand(band);
        _product.setModified(false);

        band.setData(ProductData.createInstance(new float[_sceneWidth * _sceneHeight]));
        assertEquals("data initialized, modified flag should not be set", false, _product.isModified());

        band.setData(ProductData.createInstance(new float[_sceneWidth * _sceneHeight]));
        assertEquals("data modified, modified flag should be set", true, _product.isModified());

        band.setModified(false);
        _product.setModified(false);

        band.setData(null);
        assertEquals("data set to null, modified flag should be set", true, _product.isModified());
    }

    public void testModifiedFlagDelegation() {

        Band band1 = new Band("band1", ProductData.TYPE_FLOAT32, _sceneWidth, _sceneHeight);
        Band band2 = new Band("band2", ProductData.TYPE_FLOAT32, _sceneWidth, _sceneHeight);

        _product.addBand(band1);
        _product.addBand(band2);
        _product.setModified(false);

        band1.setModified(true);
        assertEquals(true, band1.isModified());
        assertEquals(false, band2.isModified());
        assertEquals(true, _product.isModified());

        band2.setModified(true);
        assertEquals(true, band1.isModified());
        assertEquals(true, band2.isModified());
        assertEquals(true, _product.isModified());

        _product.setModified(false);
        assertEquals(false, band1.isModified());
        assertEquals(false, band2.isModified());
        assertEquals(false, _product.isModified());
    }

    public void testSetGeocoding() {
        MapProjection projection = createMapProjectionForTestSetGeocoding();
        MapInfo mapInfo = new MapInfo(projection, 0, 0, 23, 24, 12, 13, Datum.WGS_84);
        MapGeoCoding mapGeoCoding = new MapGeoCoding(mapInfo);
        int sceneRasterWidth = 243;
        int sceneRasterHeight = 524;
        Product product = new Product("name", "type", sceneRasterWidth, sceneRasterHeight);

        mapInfo.setSceneWidth(sceneRasterWidth + 1);
        mapInfo.setSceneHeight(sceneRasterHeight);
        try {
            product.setGeoCoding(mapGeoCoding);
            fail("IllegalArgumentException expected");
        } catch (IllegalArgumentException e) {
            // IllegalArgumentException expected
        }

        mapInfo.setSceneWidth(sceneRasterWidth);
        mapInfo.setSceneHeight(sceneRasterHeight + 1);

        try {
            product.setGeoCoding(mapGeoCoding);
            fail("IllegalArgumentException expected");
        } catch (IllegalArgumentException e) {
            // IllegalArgumentException expected
        }

        mapInfo.setSceneWidth(sceneRasterWidth);
        mapInfo.setSceneHeight(sceneRasterHeight);

        try {
            product.setGeoCoding(mapGeoCoding);
            // IllegalArgumentException not expected
        } catch (IllegalArgumentException e) {
            fail("IllegalArgumentException not expected");
        }
    }


    public void testUniqueGeoCodings() {
        Product p = new Product("N", "T", 4, 4);

        assertFalse(p.isUsingSingleGeoCoding());

        final GeoCoding gc1 = new ProductUtilsTest.SGeoCoding();
        final GeoCoding gc2 = new ProductUtilsTest.DGeoCoding();
        p.setGeoCoding(gc1);

        assertTrue(p.isUsingSingleGeoCoding());

        p.addBand("A", ProductData.TYPE_INT8);
        p.addBand("B", ProductData.TYPE_INT8);

        assertTrue(p.isUsingSingleGeoCoding());

        p.getBand("A").setGeoCoding(gc1);
        p.getBand("B").setGeoCoding(gc2);

        assertFalse(p.isUsingSingleGeoCoding());

        p.getBand("B").setGeoCoding(gc1);

        assertTrue(p.isUsingSingleGeoCoding());
    }

    public void testContainsPixel() {
        Product p = new Product("x", "y", 1121, 2241);

        assertTrue(p.containsPixel(0.0f, 0.0f));
        assertTrue(p.containsPixel(0.0f, 2241.0f));
        assertTrue(p.containsPixel(1121.0f, 0.0f));
        assertTrue(p.containsPixel(1121.0f, 2241.0f));
        assertTrue(p.containsPixel(500.0f, 1000.0f));

        assertFalse(p.containsPixel(-0.1f, 0.0f));
        assertFalse(p.containsPixel(0.0f, 2241.1f));
        assertFalse(p.containsPixel(1121.0f, -0.1f));
        assertFalse(p.containsPixel(1121.1f, 2241.0f));
        assertFalse(p.containsPixel(-1, -1));

        p.dispose();
    }

    public void testReadBitmaskWithByteValues() throws ParseException,
            IOException {
        Product product = new Product("Y", "X", 4, 4);
        Band band = product.addBand("flags", ProductData.TYPE_INT8);
        band.setSynthetic(true);

        final byte F1 = 0x01;
        final byte F2 = 0x02;
        final byte F3 = 0x04;

        FlagCoding flagCoding = new FlagCoding("flags");
        flagCoding.addFlag("F1", F1, null);
        flagCoding.addFlag("F2", F2, null);
        flagCoding.addFlag("F3", F3, null);

        product.addFlagCoding(flagCoding);
        band.setFlagCoding(flagCoding);

        band.ensureRasterData();
        final byte[] elems = new byte[]{
                0, F1, F2, F3,
                F1, 0, F1 + F2, F1 + F3,
                F2, F1 + F2, 0, F2 + F3,
                F3, F1 + F3, F2 + F3, 0,
        };
        band.getRasterData().setElems(elems);
        product.setModified(false);


        final byte TRUE = 23;
        final byte FALSE = 45;

        final Term termF1 = product.createTerm("flags.F1");
        final byte[] F1_MASK = new byte[]{
                FALSE, TRUE, FALSE, FALSE,
                TRUE, FALSE, TRUE, TRUE,
                FALSE, TRUE, FALSE, FALSE,
                FALSE, TRUE, FALSE, FALSE
        };
        final byte[] F1_MASK_SUB = new byte[]{
                FALSE, TRUE,
                TRUE, FALSE,
        };
        final byte[] currentF1 = new byte[4 * 4];
        final byte[] currentF1Sub = new byte[2 * 2];
        product.readBitmask(0, 0, 4, 4, termF1, currentF1, TRUE, FALSE, ProgressMonitor.NULL);
        product.readBitmask(1, 1, 2, 2, termF1, currentF1Sub, TRUE, FALSE, ProgressMonitor.NULL);
        assertEquals(true, ObjectUtils.equalObjects(currentF1, F1_MASK));
        assertEquals(true, ObjectUtils.equalObjects(currentF1Sub, F1_MASK_SUB));


        final Term termF1AF2 = product.createTerm("flags.F1 AND flags.F2");
        final byte[] F1_AND_F2_MASK = new byte[]{
                FALSE, FALSE, FALSE, FALSE,
                FALSE, FALSE, TRUE, FALSE,
                FALSE, TRUE, FALSE, FALSE,
                FALSE, FALSE, FALSE, FALSE
        };
        final byte[] F1_AND_F2_MASK_SUB = new byte[]{
                FALSE, TRUE,
                TRUE, FALSE,
        };
        final byte[] currentF1AF2 = new byte[4 * 4];
        final byte[] currentF1AF2Sub = new byte[2 * 2];
        product.readBitmask(0, 0, 4, 4, termF1AF2, currentF1AF2, TRUE, FALSE, ProgressMonitor.NULL);
        product.readBitmask(1, 1, 2, 2, termF1AF2, currentF1AF2Sub, TRUE, FALSE, ProgressMonitor.NULL);
        assertEquals(true, ObjectUtils.equalObjects(currentF1AF2, F1_AND_F2_MASK));
        assertEquals(true, ObjectUtils.equalObjects(currentF1AF2Sub, F1_AND_F2_MASK_SUB));


        final Term termF1AF2OF3 = product.createTerm("(flags.F1 AND flags.F2) OR flags.F3");
        final byte[] F1_AND_F2_OR_F3_MASK = new byte[]{
                FALSE, FALSE, FALSE, TRUE,
                FALSE, FALSE, TRUE, TRUE,
                FALSE, TRUE, FALSE, TRUE,
                TRUE, TRUE, TRUE, FALSE
        };
        final byte[] F1_AND_F2_OR_F3_MASK_SUB = new byte[]{
                FALSE, TRUE,
                TRUE, FALSE
        };
        final byte[] currentF1AF2OF3 = new byte[4 * 4];
        final byte[] currentF1AF2OF3Sub = new byte[2 * 2];
        product.readBitmask(0, 0, 4, 4, termF1AF2OF3, currentF1AF2OF3, TRUE, FALSE, ProgressMonitor.NULL);
        product.readBitmask(1, 1, 2, 2, termF1AF2OF3, currentF1AF2OF3Sub, TRUE, FALSE, ProgressMonitor.NULL);
        assertEquals(true, ObjectUtils.equalObjects(currentF1AF2OF3, F1_AND_F2_OR_F3_MASK));
        assertEquals(true, ObjectUtils.equalObjects(currentF1AF2OF3Sub, F1_AND_F2_OR_F3_MASK_SUB));
    }

    public void testReadBitmaskWithIntValues() throws ParseException,
            IOException {
        Product product = new Product("Y", "X", 4, 4);
        Band band = product.addBand("flags", ProductData.TYPE_INT8);
        band.setSynthetic(true);
        final byte F1 = 0x01;
        final byte F2 = 0x02;
        final byte F3 = 0x04;

        FlagCoding flagCoding = new FlagCoding("flags");
        flagCoding.addFlag("F1", F1, null);
        flagCoding.addFlag("F2", F2, null);
        flagCoding.addFlag("F3", F3, null);

        product.addFlagCoding(flagCoding);
        band.setFlagCoding(flagCoding);

        band.ensureRasterData();
        final byte[] elems = new byte[]{
                0, F1, F2, F3,
                F1, 0, F1 + F2, F1 + F3,
                F2, F1 + F2, 0, F2 + F3,
                F3, F1 + F3, F2 + F3, 0,
        };
        band.getRasterData().setElems(elems);
        product.setModified(false);


        final int TRUE = 23345;
        final int FALSE = 454236;

        final Term termF1 = product.createTerm("flags.F1");
        final int[] F1_MASK = new int[]{
                FALSE, TRUE, FALSE, FALSE,
                TRUE, FALSE, TRUE, TRUE,
                FALSE, TRUE, FALSE, FALSE,
                FALSE, TRUE, FALSE, FALSE
        };
        final int[] F1_MASK_SUB = new int[]{
                FALSE, TRUE,
                TRUE, FALSE,
        };
        final int[] currentF1 = new int[4 * 4];
        final int[] currentF1Sub = new int[2 * 2];
        product.readBitmask(0, 0, 4, 4, termF1, currentF1, TRUE, FALSE);
        product.readBitmask(1, 1, 2, 2, termF1, currentF1Sub, TRUE, FALSE);
        assertEquals(true, ObjectUtils.equalObjects(currentF1, F1_MASK));
        assertEquals(true, ObjectUtils.equalObjects(currentF1Sub, F1_MASK_SUB));


        final Term termF1AF2 = product.createTerm("flags.F1 AND flags.F2");
        final int[] F1_AND_F2_MASK = new int[]{
                FALSE, FALSE, FALSE, FALSE,
                FALSE, FALSE, TRUE, FALSE,
                FALSE, TRUE, FALSE, FALSE,
                FALSE, FALSE, FALSE, FALSE
        };
        final int[] F1_AND_F2_MASK_SUB = new int[]{
                FALSE, TRUE,
                TRUE, FALSE,
        };
        final int[] currentF1AF2 = new int[4 * 4];
        final int[] currentF1AF2Sub = new int[2 * 2];
        product.readBitmask(0, 0, 4, 4, termF1AF2, currentF1AF2, TRUE, FALSE);
        product.readBitmask(1, 1, 2, 2, termF1AF2, currentF1AF2Sub, TRUE, FALSE);
        assertEquals(true, ObjectUtils.equalObjects(currentF1AF2, F1_AND_F2_MASK));
        assertEquals(true, ObjectUtils.equalObjects(currentF1AF2Sub, F1_AND_F2_MASK_SUB));


        final Term termF1AF2OF3 = product.createTerm("(flags.F1 AND flags.F2) OR flags.F3");
        final int[] F1_AND_F2_OR_F3_MASK = new int[]{
                FALSE, FALSE, FALSE, TRUE,
                FALSE, FALSE, TRUE, TRUE,
                FALSE, TRUE, FALSE, TRUE,
                TRUE, TRUE, TRUE, FALSE
        };
        final int[] F1_AND_F2_OR_F3_MASK_SUB = new int[]{
                FALSE, TRUE,
                TRUE, FALSE
        };
        final int[] currentF1AF2OF3 = new int[4 * 4];
        final int[] currentF1AF2OF3Sub = new int[2 * 2];
        product.readBitmask(0, 0, 4, 4, termF1AF2OF3, currentF1AF2OF3, TRUE, FALSE);
        product.readBitmask(1, 1, 2, 2, termF1AF2OF3, currentF1AF2OF3Sub, TRUE, FALSE);
        assertEquals(true, ObjectUtils.equalObjects(currentF1AF2OF3, F1_AND_F2_OR_F3_MASK));
        assertEquals(true, ObjectUtils.equalObjects(currentF1AF2OF3Sub, F1_AND_F2_OR_F3_MASK_SUB));
    }

    public void testEnsureValidMask() throws ParseException,
            IOException {
        final Product product = new Product("n", "t", 18, 2);
        final Band flagsBand = product.addBand("flags", ProductData.TYPE_INT8);
        flagsBand.setSynthetic(true);
        final FlagCoding flagCoding = new FlagCoding("fc");
        final int f1Mask = 1;
        flagCoding.addFlag("f1", f1Mask, "");
        flagsBand.setFlagCoding(flagCoding);
        product.addFlagCoding(flagCoding);
        final byte[] elems = new byte[]{
                0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0,
                1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1,
        };
        flagsBand.setDataElems(elems);
        product.setModified(false);


        final Term term = product.createTerm("flags.f1");
        final BitRaster validMask = product.createValidMask(term, ProgressMonitor.NULL);

        for (int i = 0; i < elems.length; i++) {
            assertEquals(elems[i] == 1, validMask.isSet(i));
        }
    }

    public void testMaskProductData() throws ParseException,
            IOException {
        final Product product = new Product("Y", "X", 4, 4);
        final Band band = product.addBand("flags", ProductData.TYPE_INT8);
        band.setSynthetic(true);

        final byte F1 = 0x01;

        final FlagCoding flagCoding = new FlagCoding("flags");
        flagCoding.addFlag("F1", F1, null);

        product.addFlagCoding(flagCoding);
        band.setFlagCoding(flagCoding);

        final byte[] elems = new byte[]{
                0, F1, 0, 0,
                F1, 0, F1, F1,
                0, F1, 0, 0,
                0, F1, 0, 0,
        };
        band.setDataElems(elems);
        product.setModified(false);


        final Term term = product.createTerm("flags.F1");
        final float[] currentData = new float[4 * 4];
        final float[] currentDataSub = new float[2 * 2];
        final ProductData rasterData = ProductData.createInstance(currentData);
        final ProductData rasterDataSub = ProductData.createInstance(currentDataSub);
        final float mv = 23.345f;
        product.maskProductData(0, 0, 4, 4, term, rasterData, true, mv, ProgressMonitor.NULL);
        product.maskProductData(1, 1, 2, 2, term, rasterDataSub, true, mv, ProgressMonitor.NULL);
        final float[] expectedData = new float[]{
                0, mv, 0, 0,
                mv, 0, mv, mv,
                0, mv, 0, 0,
                0, mv, 0, 0,
        };
        final float[] expectedDataSub = new float[]{
                0, mv,
                mv, 0,
        };
        assertEquals(true, ObjectUtils.equalObjects(expectedData, currentData));
        assertEquals(true, ObjectUtils.equalObjects(expectedDataSub, currentDataSub));
    }

    public void testExpressionIsChangedIfANodeNameIsChanged() {
        final Product product = new Product("p", "t", 10, 10);
        final VirtualBand virtualBand = new VirtualBand("vb", ProductData.TYPE_FLOAT32, 10, 10,
                                                        "band1 + band2 - band3");
        product.addBand(virtualBand);
        product.addBand("band1", ProductData.TYPE_FLOAT32);
        product.addBand("band2", ProductData.TYPE_FLOAT32);
        product.addBand("band3", ProductData.TYPE_FLOAT32);

        product.getBand("band1").setName("b1");

        assertEquals("Name 'band1' is not changed",
                     "b1 + band2 - band3", virtualBand.getExpression());
    }

    public void testThatAddBandThrowExceptionIfNameIsNotUnique() {
        final Product product = new Product("p", "t", 1, 1);
        product.addBand("band1", ProductData.TYPE_FLOAT32);
        product.addTiePointGrid(new TiePointGrid("grid", 1, 1, 0, 0, 1, 1, new float[]{0.0f}));

        try {
            product.addBand("band1", ProductData.TYPE_FLOAT32);
            fail("IllegalArgumentException expected");
        } catch (IllegalArgumentException e) {
            assertTrue(e.getMessage().indexOf("name") > -1);
        }

        try {
            product.addBand("grid", ProductData.TYPE_FLOAT32);
            fail("IllegalArgumentException expected");
        } catch (IllegalArgumentException e) {
            assertTrue(e.getMessage().indexOf("name") > -1);
        }
    }

    public void testThatAddTiePointGridThrowExceptionIfNameIsNotUnique() {
        final Product product = new Product("p", "t", 1, 1);
        product.addBand("band1", ProductData.TYPE_FLOAT32);
        product.addTiePointGrid(new TiePointGrid("grid", 1, 1, 0, 0, 1, 1, new float[]{0.0f}));

        try {
            product.addTiePointGrid(new TiePointGrid("grid", 1, 1, 0, 0, 1, 1, new float[]{0.0f}));
            fail("IllegalArgumentException expected");
        } catch (IllegalArgumentException e) {
            assertTrue(e.getMessage().indexOf("name") > -1);
        }

        try {
            product.addTiePointGrid(new TiePointGrid("band1", 1, 1, 0, 0, 1, 1, new float[]{0.0f}));
            fail("IllegalArgumentException expected");
        } catch (IllegalArgumentException e) {
            assertTrue(e.getMessage().indexOf("name") > -1);
        }
    }

    public void testPreferredTileSizeProperty() {
        Product product;

        product = new Product("A", "B", 1000, 2000);
        assertEquals(null, product.getPreferredTileSize());

        product.setPreferredTileSize(new Dimension(128, 256));
        assertEquals(new Dimension(128, 256), product.getPreferredTileSize());

        product.setPreferredTileSize(new Dimension(300, 400));
        assertEquals(new Dimension(300, 400), product.getPreferredTileSize());

        product.setPreferredTileSize(null);
        assertEquals(null, product.getPreferredTileSize());
    }

    private static MapProjection createMapProjectionForTestSetGeocoding() {
        MapTransform mapTransform = new IdentityTransformDescriptor().createTransform(null);
        return new MapProjection("p1", mapTransform, "unit");
    }
}

class DummyProductReader extends AbstractProductReader {

    DummyProductReader(DummyProductReaderPlugIn plugIn) {
        super(plugIn);
    }

    @Override
    public ProductReaderPlugIn getReaderPlugIn() {
        return null;
    }

    @Override
    public Product readProductNodesImpl() throws IOException {
        throw new IllegalStateException("not implemented");
    }

    @Override
    protected void readBandRasterDataImpl(int sourceOffsetX, int sourceOffsetY, int sourceWidth, int sourceHeight,
                                          int sourceStepX, int sourceStepY, Band destBand, int destOffsetX,
                                          int destOffsetY, int destWidth, int destHeight, ProductData destBuffer,
                                          ProgressMonitor pm) throws IOException {
        throw new IllegalStateException("not implemented");
    }

    @Override
    public void close() throws IOException {
    }
}

class DummyProductReaderPlugIn implements ProductReaderPlugIn {

    public DecodeQualification getDecodeQualification(Object input) {
        return DecodeQualification.UNABLE;
    }

    public String[] getFormatNames() {
        return new String[0];
    }

    public String[] getDefaultFileExtensions() {
        return new String[0];
    }

    public Class[] getInputTypes() {
        return new Class[0];
    }

    public String getDescription(Locale locale) {
        return null;
    }

    public ProductReader createReaderInstance() {
        return new DummyProductReader(this);
    }

    public BeamFileFilter getProductFileFilter() {
        return new BeamFileFilter(getFormatNames()[0], getDefaultFileExtensions(), getDescription(null));
    }

}

class DummyProductNodeListener implements ProductNodeListener {

    public DummyProductNodeListener() {
    }

    public void nodeChanged(ProductNodeEvent event) {
    }

    public void nodeDataChanged(ProductNodeEvent event) {
    }

    public void nodeAdded(ProductNodeEvent event) {
    }

    public void nodeRemoved(ProductNodeEvent event) {
    }
}