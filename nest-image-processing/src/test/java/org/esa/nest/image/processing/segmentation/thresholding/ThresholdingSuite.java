/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.esa.nest.image.processing.segmentation.thresholding;

import org.esa.nest.image.processing.segmentation.thresholding.BasicThresholdingOpTest;
import org.esa.nest.image.processing.segmentation.thresholding.BasicThresholdingOpUITest;
import org.esa.nest.image.processing.segmentation.thresholding.ThresholdingMethodEnforcementTest;
import org.esa.nest.image.processing.segmentation.thresholding.ThresholdingTypeOperatorTest;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.runner.RunWith;
import org.junit.runners.Suite;

/**
 *
 * @author Emanuela
 */
@RunWith(Suite.class)
@Suite.SuiteClasses({ThresholdingTypeOperatorTest.class,
    ThresholdingMethodEnforcementTest.class,
    BasicThresholdingOpTest.class,
    BasicThresholdingOpUITest.class})
public class ThresholdingSuite {

    @BeforeClass
    public static void setUpClass() throws Exception {
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
    }

    @Before
    public void setUp() throws Exception {
    }

    @After
    public void tearDown() throws Exception {
    }
}
