/*
 * $Id: GETASSE30ElevationModelDescriptorTest.java,v 1.1 2009-04-28 14:37:14 lveci Exp $
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
package org.esa.beam.dataio.getasse30;

import junit.framework.TestCase;

public class GETASSE30ElevationModelDescriptorTest extends TestCase {

    final GETASSE30ElevationModelDescriptor _descriptor = new GETASSE30ElevationModelDescriptor();

    public void testConstantProperties() {
        assertEquals("GETASSE30", _descriptor.getName());
    }

    public void testFilenameCreation() {

        assertEquals("45S004W.GETASSE30", _descriptor.createTileFilename(-45, -4));
        assertEquals("45S004E.GETASSE30", _descriptor.createTileFilename(-45, +4));
        assertEquals("45N004W.GETASSE30", _descriptor.createTileFilename(+45, -4));
        assertEquals("45N004E.GETASSE30", _descriptor.createTileFilename(+45, +4));

        assertEquals("05S045W.GETASSE30", _descriptor.createTileFilename(-5, -45));
        assertEquals("05S045E.GETASSE30", _descriptor.createTileFilename(-5, +45));
        assertEquals("05N045W.GETASSE30", _descriptor.createTileFilename(+5, -45));
        assertEquals("05N045E.GETASSE30", _descriptor.createTileFilename(+5, +45));

        assertEquals("90S180W.GETASSE30", _descriptor.createTileFilename(-90, -180));
        assertEquals("90S180E.GETASSE30", _descriptor.createTileFilename(-90, +180));
        assertEquals("90N180W.GETASSE30", _descriptor.createTileFilename(+90, -180));
        assertEquals("90N180E.GETASSE30", _descriptor.createTileFilename(+90, +180));
    }

}
