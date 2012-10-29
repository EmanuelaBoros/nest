/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.esa.nest.image.processing.segmentation.thresholding;

import com.bc.ceres.core.ProgressMonitor;
import ij.process.ByteProcessor;
import java.util.Map;
import org.esa.beam.framework.datamodel.Band;
import org.esa.beam.framework.gpf.Tile;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import static org.junit.Assert.*;

/**
 *
 * @author Emanuela
 */
public class ThresholdingTypeOperatorTest {
    
    public ThresholdingTypeOperatorTest() {
    }
    
    @BeforeClass
    public static void setUpClass() {
    }
    
    @AfterClass
    public static void tearDownClass() {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }

    /**
     * Test of values method, of class ThresholdingTypeOperator.
     */
//    @Test
    public void testValues() {
        System.out.println("values");
        ThresholdingTypeOperator[] expResult = null;
        ThresholdingTypeOperator[] result = ThresholdingTypeOperator.values();
        assertArrayEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of valueOf method, of class ThresholdingTypeOperator.
     */
//    @Test
    public void testValueOf() {
        System.out.println("valueOf");
        String string = "";
        ThresholdingTypeOperator expResult = null;
        ThresholdingTypeOperator result = ThresholdingTypeOperator.valueOf(string);
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of computeThresholdingOperator method, of class ThresholdingTypeOperator.
     */
//    @Test
    public void testComputeThresholdingOperator() {
        System.out.println("computeThresholdingOperator");
        Band band = null;
        Tile tile = null;
        Tile tile1 = null;
        int i = 0;
        int i1 = 0;
        int i2 = 0;
        int i3 = 0;
        ProgressMonitor pm = null;
        Map<String, Object> map = null;
        ThresholdingTypeOperator instance = null;
        ByteProcessor expResult = null;
        ByteProcessor result = instance.computeThresholdingOperator(band, tile, tile1, i, i1, i2, i3, pm, map);
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

//    public class ThresholdingTypeOperatorImpl extends ThresholdingTypeOperator {
//
//        public ByteProcessor computeThresholdingOperator(Band band, Tile tile, Tile tile1, int i, int i1, int i2, int i3, ProgressMonitor pm, Map<String, Object> map) {
//            return null;
//        }
//    }
}
