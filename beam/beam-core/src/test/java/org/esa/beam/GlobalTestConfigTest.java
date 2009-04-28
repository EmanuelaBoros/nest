/*
 * $Id: GlobalTestConfigTest.java,v 1.1 2009-04-28 14:39:33 lveci Exp $
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
package org.esa.beam;

import java.io.File;
import java.util.Properties;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

import org.esa.beam.util.SystemUtils;

public class GlobalTestConfigTest extends TestCase {

    private Properties _propertys;

    public GlobalTestConfigTest(String testName) {
        super(testName);
    }

    public static Test suite() {
        return new TestSuite(GlobalTestConfigTest.class);
    }

    @Override
    protected void setUp() throws Exception {
        _propertys = System.getProperties();
    }

    @Override
    protected void tearDown() throws Exception {
        System.setProperties(_propertys);
    }

    public void testGetEnvisatTestDataDir() {
        File file = null;
        try {
            file = GlobalTestConfig.getBeamTestDataInputDirectory();
        } catch (SecurityException e) {
            fail("SecurityException not expected");
        }
        assertNotNull(file);

        System.setProperty(GlobalTestConfig.BEAM_TEST_DATA_INPUT_DIR_PROPERTY_NAME,
                           SystemUtils.convertToLocalPath("C:/envi/test/data/"));
        try {
            file = GlobalTestConfig.getBeamTestDataInputDirectory();
        } catch (SecurityException e) {
            fail("SecurityException not expected");
        }
        assertEquals(new File(SystemUtils.convertToLocalPath("C:/envi/test/data/")), file);

        System.getProperties().remove(GlobalTestConfig.BEAM_TEST_DATA_INPUT_DIR_PROPERTY_NAME);
        try {
            file = GlobalTestConfig.getBeamTestDataInputDirectory();
        } catch (SecurityException e) {
            fail("SecurityException not expected");
        }
        final File defaultFile = new File(SystemUtils.getBeamHomeDir(),
                                          SystemUtils.convertToLocalPath(
                                                  GlobalTestConfig.BEAM_TEST_DATA_INPUT_DIR_DEFAULT_PATH));
        assertEquals(defaultFile, file);
    }

}
