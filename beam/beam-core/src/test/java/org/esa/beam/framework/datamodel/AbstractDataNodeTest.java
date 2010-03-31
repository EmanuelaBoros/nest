/*
 * $Id: AbstractDataNodeTest.java,v 1.2 2010-03-31 13:59:56 lveci Exp $
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


public class AbstractDataNodeTest extends AbstractNamedNodeTest {

    public AbstractDataNodeTest(String testName) {
        super(testName);
    }    

    @Override
    public void testSetUnit(DataNode dataNode) {

        // old value --> null ?
        dataNode.setUnit(null);
        assertEquals(null, dataNode.getUnit());
        assertEquals(false, dataNode.isModified());

        // null --> new value: is modified ?
        dataNode.setUnit("mg/m^3");
        assertEquals("mg/m^3", dataNode.getUnit());
        assertEquals(true, dataNode.isModified());

        // old value == new value?
        dataNode.setModified(false);
        dataNode.setUnit("mg/m^3");
        assertEquals("mg/m^3", dataNode.getUnit());
        assertEquals(false, dataNode.isModified());

        // old value != new value?
        dataNode.setUnit("g/cm^3");
        assertEquals("g/cm^3", dataNode.getUnit());
        assertEquals(true, dataNode.isModified());
    }
}
