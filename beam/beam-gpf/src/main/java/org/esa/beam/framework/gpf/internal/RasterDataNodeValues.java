/*
 * $Id: RasterDataNodeValues.java,v 1.1 2009-12-02 16:52:11 lveci Exp $
 * 
 * Copyright (C) 2009 by Brockmann Consult (info@brockmann-consult.de)
 * 
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation. This program is distributed in the hope it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 * 
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place - Suite 330, Boston, MA 02111-1307, USA.
 */
package org.esa.beam.framework.gpf.internal;

import org.esa.beam.framework.datamodel.Band;
import org.esa.beam.framework.datamodel.Mask;
import org.esa.beam.framework.datamodel.Product;
import org.esa.beam.framework.datamodel.RasterDataNode;
import org.esa.beam.framework.datamodel.TiePointGrid;

/**
 * Returns the names of all elements of the given product that are of the given
 * subtype of {@code ProductNode}.
 * 
 * @author Marco Zuehlke
 * @since BEAM 4.7
 */
public class RasterDataNodeValues {
    
    public static final String ATTRIBUTE_NAME = "rasterDataNodeType";

    /**
     * Returns the names of all elements of the given product that are of the
     * given subtype of {@code RasterDataNode}.
     * 
     * @param product
     *            The product.
     * @param rasterDataNodeType
     *            The subtype of {@code RasterDataNode}.
     * @return the names of all matching elements
     */
    public static String[] getNames(Product product, Class<? extends RasterDataNode> rasterDataNodeType) {
        return getNames(product, rasterDataNodeType, false);
    }

    /**
     * Returns the names of all elements of the given product that are of the
     * given subtype of {@code RasterDataNode}. Optionally an empty value is included.
     * 
     * @param product
     *            The product.
     * @param rasterDataNodeType
     *            The subtype of {@code RasterDataNode}.
     * @param includeEmptyValue
     *            If {code true} and empty {@code String} is included as the first value.
     * @return the names of all matching elements
     */
    public static String[] getNames(Product product, Class<? extends RasterDataNode> rasterDataNodeType,
                                          boolean includeEmptyValue) {
        String[] allValues;
        if (includeEmptyValue) {
            String[] valueNames = getValueNames(product, rasterDataNodeType);
            allValues = new String[valueNames.length + 1];
            allValues[0] = "";
            System.arraycopy(valueNames, 0, allValues, 1, valueNames.length);
        } else {
            allValues = getValueNames(product, rasterDataNodeType);
        }
        return allValues;
    }

    private static String[] getValueNames(Product product, Class<? extends RasterDataNode> rasterDataNodeType) {
        if (rasterDataNodeType == Band.class) {
            return product.getBandNames();
        } else if (rasterDataNodeType == Mask.class) {
            return product.getMaskGroup().getNodeNames();
        } else if (rasterDataNodeType == TiePointGrid.class) {
            return product.getTiePointGridNames();            
        }
        throw new IllegalArgumentException("Unsupported 'rasterDataNodeType': " + rasterDataNodeType);
    }
}
