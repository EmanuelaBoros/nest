/*
 * $Id: ProductIOTest.java,v 1.1 2009-04-28 14:39:33 lveci Exp $
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

package org.esa.beam.framework.dataio;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

import java.io.File;
import java.io.IOException;

public class ProductIOTest extends TestCase {

    public ProductIOTest(String testName) {
        super(testName);
    }

    public static Test suite() {
        return new TestSuite(ProductIOTest.class);
    }

    public void testThatDefaultReaderAndWriterAreImplemented() {
        assertNotNull(ProductIO.getProductReader("BEAM-DIMAP"));
        assertNotNull(ProductIO.getProductWriter("BEAM-DIMAP"));
    }

    public void testReadProductArgsChecking() {
        try {
            ProductIO.readProduct((File) null, null);
            fail();
        } catch (IOException expected) {
            fail();
        } catch (IllegalArgumentException expected) {
        }

        try {
            ProductIO.readProduct("rallala", null);
            fail();
        } catch (IOException expected) {
        }
    }

}