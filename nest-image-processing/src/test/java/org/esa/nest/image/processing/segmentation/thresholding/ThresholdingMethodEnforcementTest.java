/*
 * Copyright (C) 2012 by Array Systems Computing Inc. http://www.array.ca
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation; either version 3 of the License, or (at your option)
 * any later version.
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, see http://www.gnu.org/licenses/
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
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author Emanuela
 */
public class ThresholdingMethodEnforcementTest {

    public ThresholdingMethodEnforcementTest() {
    }

    @Before
    public void setUp() {
    }

    @After
    public void tearDown() {
    }

    /**
     * Test of computeThresholdingOperator method, of class
     * ThresholdingMethodEnforcement.
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
        ThresholdingMethodEnforcement instance = new ThresholdingMethodEnforcementImpl();
        ByteProcessor expResult = null;
        ByteProcessor result = instance.computeThresholdingOperator(band, tile, tile1, i, i1, i2, i3, pm, map);
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of getThresholdedImage method, of class
     * ThresholdingMethodEnforcement.
     */
//    @Test
    public void testGetThresholdedImage() {
        System.out.println("getThresholdedImage");
        ByteProcessor bp = null;
        ThresholdingMethodEnforcement instance = new ThresholdingMethodEnforcementImpl();
        ByteProcessor expResult = null;
        ByteProcessor result = instance.getThresholdedImage(bp);
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    public class ThresholdingMethodEnforcementImpl implements ThresholdingMethodEnforcement {

        public ByteProcessor computeThresholdingOperator(Band band, Tile tile, Tile tile1, int i, int i1, int i2, int i3, ProgressMonitor pm, Map<String, Object> map) {
            return null;
        }

        public ByteProcessor getThresholdedImage(ByteProcessor bp) {
            return null;
        }
    }
}
