/*
 * $Id: ProductNodeFilter.java,v 1.1 2009-04-28 14:39:33 lveci Exp $
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
 * A filter for abstract product nodes.
 * <p/>
 * <p> Instances of this interface may be passed to the {@link ProductNodeList#createSubset(ProductNodeFilter)} method.
 */
public interface ProductNodeFilter {

    /**
     * Tests whether or not the specified abstract product node should be included in a product node list.
     *
     * @param productNode the product node to be tested
     *
     * @return true if and only if product node should be included
     */
    boolean accept(ProductNode productNode);
}
