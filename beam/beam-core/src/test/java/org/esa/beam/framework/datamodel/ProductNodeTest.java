/*
 * $Id: ProductNodeTest.java,v 1.1 2009-04-28 14:39:33 lveci Exp $
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

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

public class ProductNodeTest extends TestCase {

    public ProductNodeTest(String testName) {
        super(testName);
    }

    public static Test suite() {
        return new TestSuite(ProductNodeTest.class);
    }

    @Override
    protected void setUp() {
    }

    public void testSetOwnerToNullAfterNodeRemoval() {
        final ProductNode[] owner = new Product[1];
        final Product product = new Product("product", "t", 10, 10);
        product.addProductNodeListener(new ProductNodeListenerAdapter() {
            @Override
            public void nodeRemoved(ProductNodeEvent event) {
                ProductNode sourceNode = event.getSourceNode();
                owner[0] = sourceNode.getOwner();
            }
        });

        product.addBand("band1", ProductData.TYPE_INT16);
        Band addedBand = product.getBandAt(0);
        assertSame(product, addedBand.getOwner());

        product.removeBand(addedBand);
        // ensures that the given source node
        // in the node removed event has an owner
        // at the time the event fired.
        assertSame(product, owner[0]);
        // After all listeners notified, that the
        // node was removed, the owner of the node
        // must be set to null
        assertNull(addedBand.getOwner());
    }

    public void testSetProductToNullAfterNodeRemoval() {

        Band band = new Band("b", ProductData.TYPE_INT16, 10, 10);

        assertNull(band.getOwner());
        assertNull(band.getProduct());

        final Product p1 = new Product("p1", "t", 10, 10);
        p1.addBand(band);

        assertSame(p1, band.getOwner());
        assertSame(p1, band.getProduct());

        p1.removeBand(band);

        assertNull(band.getOwner());
        assertNull(band.getProduct());

        final Product p2 = new Product("p2", "t", 10, 10);
        p2.addBand(band);

        assertSame(p2, band.getOwner());
        assertSame(p2, band.getProduct());

        p2.removeBand(band);

        assertNull(band.getOwner());
        assertNull(band.getProduct());
    }

    public void testSetName() {
        int numberOfExceptionsTrown;
        int expectedNumberOfExceptions;
        final ProductNode productNode = new Band("valid", ProductData.TYPE_INT8, 1, 1);


        final String[] invalidNames = new String[]{
            ".Band", "",  " ",  "       ", "or", "not", "and"
        };
        numberOfExceptionsTrown = tryToSetInvalidNames(productNode, invalidNames);
        expectedNumberOfExceptions = invalidNames.length;
        assertEquals(expectedNumberOfExceptions, numberOfExceptionsTrown);


        final String[] validNames = new String[]{
            "Band", "band1", "ba_nd", "band_", "_band", "Band.sdf",
            "1band", "ba#nd", "band~", "band 2", "band ", " band"
        };
        numberOfExceptionsTrown = tryToSetValidNames(productNode, validNames);
        expectedNumberOfExceptions = 0;
        assertEquals(expectedNumberOfExceptions, numberOfExceptionsTrown);
    }

    public void testIsValidNodeName(){
        assertFalse(ProductNode.isValidNodeName(""));
        assertFalse(ProductNode.isValidNodeName(" "));
        assertFalse(ProductNode.isValidNodeName("\\"));
        assertFalse(ProductNode.isValidNodeName("/"));
        assertFalse(ProductNode.isValidNodeName("*"));
        assertFalse(ProductNode.isValidNodeName("?"));
        assertFalse(ProductNode.isValidNodeName("\""));
        assertFalse(ProductNode.isValidNodeName(":"));
        assertFalse(ProductNode.isValidNodeName("<"));
        assertFalse(ProductNode.isValidNodeName(">"));
        assertFalse(ProductNode.isValidNodeName("|"));
        assertFalse(ProductNode.isValidNodeName("."));
        assertFalse(ProductNode.isValidNodeName(".a"));

        assertTrue(ProductNode.isValidNodeName("a"));
        assertTrue(ProductNode.isValidNodeName("a."));
        assertTrue(ProductNode.isValidNodeName("1"));
        assertTrue(ProductNode.isValidNodeName("F"));
        assertTrue(ProductNode.isValidNodeName("_"));
        assertTrue(ProductNode.isValidNodeName("-"));
        assertTrue(ProductNode.isValidNodeName("$"));
        assertTrue(ProductNode.isValidNodeName("\u20ac")); // Euro
        assertTrue(ProductNode.isValidNodeName("@"));
        assertTrue(ProductNode.isValidNodeName("+"));
        assertTrue(ProductNode.isValidNodeName("~"));
    }

    private int tryToSetInvalidNames(final ProductNode productNode, final String[] names) {
        int countedExceptions = 0;
        for (int i = 0; i < names.length; i++) {
            try {
                productNode.setName(names[i]);
                fail("IllegalArgumentException expected for name '" + names[i] + "'");
            } catch (IllegalArgumentException e) {
                countedExceptions++;
            }
        }
        return countedExceptions;
    }

    private int tryToSetValidNames(final ProductNode productNode, final String[] names) {
        int countedExceptions = 0;
        for (int i = 0; i < names.length; i++) {
            try {
                productNode.setName(names[i]);
            } catch (IllegalArgumentException e) {
                countedExceptions++;
                fail("IllegalArgumentException was NOT expected for name '" + names[i] + "'");
            }
        }
        return countedExceptions;
    }
}

