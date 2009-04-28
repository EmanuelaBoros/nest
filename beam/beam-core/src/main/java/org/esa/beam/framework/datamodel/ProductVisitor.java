/*
 * $Id: ProductVisitor.java,v 1.1 2009-04-28 14:39:33 lveci Exp $
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


/**
 * A visitor for a product and all other product nodes. This interface is part of the <i>visitor pattern</i> used to
 * visit all nodes of a data product. Implementations of this interface can be passed to the <code>acceptVisitor</code>
 * method of an <code>Product</code> (or any other <code>ProductNode</code>).
 *
 * @author Norman Fomferra
 * @version $Revision: 1.1 $ $Date: 2009-04-28 14:39:33 $
 * @see org.esa.beam.framework.datamodel.Product#acceptVisitor(ProductVisitor)
 * @see ProductNode#acceptVisitor(ProductVisitor)
 */
public interface ProductVisitor {

    /**
     * Visits a product.
     *
     * @param product the product to be visited
     */
    void visit(Product product);

    /**
     * Visits a group whithin a product.
     *
     * @param group the group to be visited
     */
    void visit(MetadataElement group);

    /**
     * Visits a tie-point grid within a product or group.
     *
     * @param grid the tie-point grid to be visited
     */
    void visit(TiePointGrid grid);

    /**
     * Visits a band within a product or group.
     *
     * @param band the band to be visited
     */
    void visit(Band band);

    /**
     * Visits a virtual band.
     *
     * @param virtualBand the bitmask definition to be visited
     */
    void visit(VirtualBand virtualBand);

    /**
     * Visits an attribute.
     *
     * @param attribute the attribute to be visited
     */
    void visit(MetadataAttribute attribute);

    /**
     * Visits a flag coding.
     *
     * @param flagCoding the flag coding to be visited
     */
    void visit(FlagCoding flagCoding);

    /**
     * Visits an index coding.
     *
     * @param indexCoding the index coding to be visited
     */
    void visit(IndexCoding indexCoding);

    /**
     * Visits a bitmask definition.
     *
     * @param bitmaskDef the bitmask definition to be visited
     */
    void visit(BitmaskDef bitmaskDef);

    /**
     * Visits a node group.
     *
     * @param group the group to be visited
     */
    void visit(ProductNodeGroup group);
}
