/*
 * $Id: ModisTiePointDescriptionTest.java,v 1.2 2010-03-31 13:59:56 lveci Exp $
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
 *//*
 * $Id: ModisTiePointDescriptionTest.java,v 1.2 2010-03-31 13:59:56 lveci Exp $
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
package org.esa.beam.dataio.modis.productdb;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

public class ModisTiePointDescriptionTest extends TestCase {

    public ModisTiePointDescriptionTest(String name) {
        super(name);
    }

    public static Test suite() {
        return new TestSuite(ModisTiePointDescriptionTest.class);
    }

    public void testTheFunctionality() {
        String expName = "tie_point_name";
        String expScale = "scale_name";
        String expOffset = "offset_name";
        String expUnit = "unit_name";

        // all value null allowed
        ModisTiePointDescription desc = new ModisTiePointDescription(null, null, null, null);
        assertEquals(null, desc.getName());
        assertEquals(null, desc.getScaleAttribName());
        assertEquals(null, desc.getOffsetAttribName());
        assertEquals(null, desc.getUnitAttribName());

        // check values one after the other
        desc = new ModisTiePointDescription(expName, null, null, null);
        assertEquals(expName, desc.getName());
        assertEquals(null, desc.getScaleAttribName());
        assertEquals(null, desc.getOffsetAttribName());
        assertEquals(null, desc.getUnitAttribName());

        desc = new ModisTiePointDescription(null, expScale, null, null);
        assertEquals(null, desc.getName());
        assertEquals(expScale, desc.getScaleAttribName());
        assertEquals(null, desc.getOffsetAttribName());
        assertEquals(null, desc.getUnitAttribName());

        desc = new ModisTiePointDescription(null, null, expOffset, null);
        assertEquals(null, desc.getName());
        assertEquals(null, desc.getScaleAttribName());
        assertEquals(expOffset, desc.getOffsetAttribName());
        assertEquals(null, desc.getUnitAttribName());

        desc = new ModisTiePointDescription(null, null, null, expUnit);
        assertEquals(null, desc.getName());
        assertEquals(null, desc.getScaleAttribName());
        assertEquals(null, desc.getOffsetAttribName());
        assertEquals(expUnit, desc.getUnitAttribName());
    }
}
