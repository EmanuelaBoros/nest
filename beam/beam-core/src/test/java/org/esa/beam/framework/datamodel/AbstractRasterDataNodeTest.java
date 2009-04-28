/*
 * $Id: AbstractRasterDataNodeTest.java,v 1.1 2009-04-28 14:39:33 lveci Exp $
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


public abstract class AbstractRasterDataNodeTest extends AbstractDataNodeTest {

    public AbstractRasterDataNodeTest(String testName) {
        super(testName);
    }

    protected abstract RasterDataNode createRasterDataNode();

    public void testSetAndGetBandStatistics() {
        RasterDataNode rasterDataNode = createRasterDataNode();
        assertEquals(null, rasterDataNode.getImageInfo());
        final ImageInfo imageInfo = new ImageInfo(new ColorPaletteDef(0, 1));
        rasterDataNode.setImageInfo(imageInfo);
        assertSame(imageInfo, rasterDataNode.getImageInfo());
    }

    public void testValidMaskExpressionIsAdjustedIfNodeNameChanged() {
        final RasterDataNode rasterDataNode = createRasterDataNode();
        rasterDataNode.setValidPixelExpression("flagsBand.f1 || not flagsBand.f2");
        final Product product = new Product("p", "NoType",
                                            rasterDataNode.getRasterWidth(), rasterDataNode.getRasterHeight());
        addRasterDataNodeToProduct(product, rasterDataNode);

        final FlagCoding flagCoding = new FlagCoding("f");
        flagCoding.addFlag("f1", 0x01, "descr");
        flagCoding.addFlag("f2", 0x02, "descr");

        final Band flagsBand = product.addBand("flagsBand", ProductData.TYPE_INT8);
        flagsBand.setFlagCoding(flagCoding);
        product.addFlagCoding(flagCoding);

        flagsBand.setName("flags");

        final String currentExpression = rasterDataNode.getValidPixelExpression();
        final String expectedExpression = "flags.f1 || not flags.f2";
        assertEquals("name is not changed", expectedExpression, currentExpression);
    }

    public void testUpdateExpression() {
        final String oldIdentifier = "oldIdent";
        final String newIdentifier = "newIdent";
        final String initialExpression = "ident_1 + oldIdent - ident_3";
        final String renamedExpression = "ident_1 + newIdent - ident_3";
        final RasterDataNode node = createRasterDataNode();
        node.setValidPixelExpression(initialExpression);
        final int width = node.getSceneRasterWidth();
        final int height = node.getSceneRasterHeight();
        final boolean[] isActiv = {false};
        final Product product = new Product("n", "t", width, height) {
            @Override
            protected void fireNodeAdded(ProductNode sourceNode) {
                if (isActiv[0]) {
                    fail("Event not expected.");
                }
            }

            @Override
            protected void fireNodeChanged(ProductNode sourceNode, String propertyName, Object oldValue) {
                if (isActiv[0] && !ProductNode.PROPERTY_NAME_MODIFIED.equalsIgnoreCase(propertyName)
                && !RasterDataNode.PROPERTY_NAME_ROI_DEFINITION.equalsIgnoreCase(propertyName)) {
                    fail("Event for property '" + propertyName + "' not expected.");
                }
            }

            @Override
            protected void fireNodeDataChanged(DataNode sourceNode) {
                if (isActiv[0]) {
                    fail("Event not expected.");
                }
            }

            @Override
            protected void fireNodeRemoved(ProductNode sourceNode) {
                if (isActiv[0]) {
                    fail("Event not expected.");
                }
            }
        };
        addRasterDataNodeToProduct(product, node);
        product.setModified(false);
        assertFalse(node.isModified());
        assertEquals(initialExpression, node.getValidPixelExpression());

        isActiv[0] = true;
        node.updateExpression(oldIdentifier, newIdentifier);
        assertTrue(node.isModified());
        assertEquals(renamedExpression, node.getValidPixelExpression());
    }

    private void addRasterDataNodeToProduct(final Product product, final RasterDataNode rasterDataNode) {
        if (rasterDataNode instanceof Band) {
            product.addBand((Band) rasterDataNode);
        } else if (rasterDataNode instanceof TiePointGrid) {
            product.addTiePointGrid((TiePointGrid) rasterDataNode);
        } else {
            fail("couldn't add RasterDataNode to product. Node is of unknown type.");
        }
    }
}
