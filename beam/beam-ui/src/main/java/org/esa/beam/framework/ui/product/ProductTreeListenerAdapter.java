/*
 * $Id: ProductTreeListenerAdapter.java,v 1.2 2009-12-11 20:46:14 lveci Exp $
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
package org.esa.beam.framework.ui.product;

import org.esa.beam.framework.datamodel.Band;
import org.esa.beam.framework.datamodel.MetadataElement;
import org.esa.beam.framework.datamodel.Product;
import org.esa.beam.framework.datamodel.ProductNode;
import org.esa.beam.framework.datamodel.TiePointGrid;
import org.esa.beam.framework.datamodel.VectorDataNode;

// todo - this is a stupid interface, use SelectionService/SelectionProvide instead (nf, 10.2009)


/**
 * A listener adapter for product tree changes.
 *
 * @since BEAM 4.7
 */
public abstract class ProductTreeListenerAdapter implements ProductTreeListener2 {
    /**
     * Called when a product has been added to the tree.
     *
     * @param product The product.
     */
    @Override
    public void productAdded(Product product) {
    }

    /**
     * Called when a product has been removed from the tree.
     *
     * @param product The product.
     */
    @Override
    public void productRemoved(Product product) {
    }

    /**
     * Called when a product has been selected in the tree.
     *
     * @param product    The selected product.
     * @param clickCount The number of mouse clicks.
     */
    @Override
    public void productSelected(Product product, int clickCount) {
    }

    /**
     * Called when a product's metadata element has been selected in the tree.
     *
     * @param matadataElement The selected metadata element.
     * @param clickCount      The number of mouse clicks.
     */
    @Override
    public void metadataElementSelected(MetadataElement matadataElement, int clickCount) {
    }

    /**
     * Called when a product's tie-point grid has been selected in the tree.
     *
     * @param tiePointGrid The selected tie-point grid.
     * @param clickCount   The number of mouse clicks.
     */
    @Override
    public void tiePointGridSelected(TiePointGrid tiePointGrid, int clickCount) {
    }

    /**
     * Called when a product's band has been selected in the tree.
     *
     * @param band       The selected band.
     * @param clickCount The number of mouse clicks.
     */
    @Override
    public void bandSelected(Band band, int clickCount) {
    }

    /**
     * Called when a product's node has been selected in the tree.
     *
     * @param vectorDataNode The selected vector data.
     * @param clickCount The number of mouse clicks.
     *
     * @since BEAM 4.7
     */
    @Override
    public void vectorDataSelected(VectorDataNode vectorDataNode, int clickCount) {
    }

    /**
     * Called when a product's node has been selected in the tree.
     *
     * @param productNode The selected product node.
     * @param clickCount  The number of mouse clicks.
     *
     * @since BEAM 4.7
     */
    @Override
    public void productNodeSelected(ProductNode productNode, int clickCount) {
    }
}