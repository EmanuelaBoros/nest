/*
 * $Id: DataItemInfo.java,v 1.1 2009-04-28 14:37:13 lveci Exp $
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
package org.esa.beam.dataio.envisat;


import org.esa.beam.framework.datamodel.ProductData;
import org.esa.beam.util.Debug;

/**
 * The <code>DataItemInfo</code> class represents a named item having a specific data type.
 *
 * @author Norman Fomferra
 * @version $Revision: 1.1 $ $Date: 2009-04-28 14:37:13 $
 */
public abstract class DataItemInfo extends ItemInfo {

    /**
     * This items's data type. Always one of the multiple <code>org.esa.beam.framework.datamodel.ProductData.TYPE_</code>X
     * constants
     */
    private final int _dataType;

    /**
     * This items's unit (optional).
     */
    private final String _physicalUnit;

    /**
     * Constructs a new data item info from the supplied parameters.
     *
     * @param itemName     the item name, must not be null or empty
     * @param dataType     the internal data type. Must be one of the multiple <code>org.esa.beam.framework.datamodel.ProductData.TYPE_</code>X
     *                     constants
     * @param physicalUnit the item's physical unit (optional, can be null)
     * @param description  the item's description (optional, can be null)
     *
     * @see org.esa.beam.framework.datamodel.ProductData
     */
    protected DataItemInfo(String itemName,
                           int dataType,
                           String physicalUnit,
                           String description) {
        super(itemName, description);
        Debug.assertTrue(dataType != ProductData.TYPE_UNDEFINED,
                         "undefined field data type"); /*I18N*/
        _dataType = dataType;
        _physicalUnit = physicalUnit;
    }

    /**
     * Gets the field's internal data type which is always one of the multiple <code>TYPE_</code>X constants defined in
     * the <code>org.esa.beam.framework.datamodel.ProductData</code> interface.
     *
     * @return the data type
     *
     * @see org.esa.beam.framework.datamodel.ProductData
     */
    public final int getDataType() {
        return _dataType;
    }

    /**
     * Gets the physical units string, which can be <code>null</code>.
     *
     * @return the physical unit
     */
    public final String getPhysicalUnit() {
        return _physicalUnit;
    }


    /**
     * Utility method which returns the size in bytes required to store a single element of the type given by the
     * supplied data type ID. If the given type is unknown, the method returns zero.
     * <p/>
     * <p>IMPORTANT NOTE: This method returns <code>12</code> (= 3 x 4 bytes) for the data type
     * <code>ProductData.TYPE_UTC</code>, since the DDDB interprets an UTC value as a single element, where as the
     * <code>ProductData.UTC</code> stores it as three <code>int</code>s.
     *
     * @param itemDataType the item's data type, must be one of the <code>org.esa.beam.framework.datamodel.ProductData.TYPE_</code>X
     *                     constants
     *
     * @return the data element size in bytes
     *
     * @see org.esa.beam.framework.datamodel.ProductData
     * @see org.esa.beam.framework.datamodel.ProductData.UTC
     */
    public static int getDataTypeElemSize(int itemDataType) {
        switch (itemDataType) {
        case ProductData.TYPE_INT8:
        case ProductData.TYPE_ASCII:
        case ProductData.TYPE_UINT8:
            return 1;
        case ProductData.TYPE_INT16:
        case ProductData.TYPE_UINT16:
            return 2;
        case ProductData.TYPE_INT32:
        case ProductData.TYPE_UINT32:
        case ProductData.TYPE_FLOAT32:
            return 4;
        case ProductData.TYPE_FLOAT64:
            return 8;
        case ProductData.TYPE_UTC:
            return 3 * 4;
        default:
            return 0;
        }
    }

}
