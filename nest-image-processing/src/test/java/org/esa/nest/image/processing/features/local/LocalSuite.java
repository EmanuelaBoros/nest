/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.esa.nest.image.processing.features.local;

import org.esa.nest.image.processing.features.local.matching.MatchingSuite;
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
@Suite.SuiteClasses({SIFTKeypointOpTest.class, MatchingSuite.class,
    ASIFTKeypointOpTest.class,
    ColourSIFTKeypointOpTest.class})
public class LocalSuite {

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
