/*
 * $Id: FieldInfo.java,v 1.1 2009-04-28 14:37:13 lveci Exp $
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
 * The <code>FieldInfo</code> class contains the information about the structure of a particular record field.
 * <p/>
 * <p> A <code>RecordInfo</code> instance contains a list of <code>FieldInfo</code> instances describing each field of
 * the record.
 *
 * @author Norman Fomferra
 * @version $Revision: 1.1 $ $Date: 2009-04-28 14:37:13 $
 * @see org.esa.beam.dataio.envisat.Field
 * @see org.esa.beam.dataio.envisat.Record
 * @see org.esa.beam.dataio.envisat.RecordInfo
 */
public class FieldInfo extends DataItemInfo {

    /**
     * The number of data elements contained in this field (field-width).
     */
    private final int _numDataElems;

    /**
     * Constructs a new field-info from the supplied parameters.
     *
     * @param fieldName    the field name, must not be null or empty
     * @param dataType     the internal data type. Must be one of the multiple <code>TYPE_</code>XXX constants defined
     *                     in the <code>org.esa.beam.framework.datamodel.ProductData</code> interface.
     * @param numDataElems the number of data elements contained in this field (field-width), must be <code>&gt;=
     *                     1</code>.
     * @param physicalUnit the field's unit (optional, can be null)
     * @param description  the field's description (optional, can be null)
     *
     * @see org.esa.beam.framework.datamodel.ProductData
     */
    FieldInfo(String fieldName,
              int dataType,
              int numDataElems,
              String physicalUnit,
              String description) {
        super(fieldName, dataType, physicalUnit, description);
        Debug.assertTrue(numDataElems >= 1,
                         "number of data elements must be greater zero"); /*I18N*/
        _numDataElems = numDataElems;
    }

    /**
     * Factory method which creates a new field instance according to the field structure defined by this field-info.
     * <p> The method simply calls <code>Field.create(this)</code> to create the instance.
     *
     * @return a new field instance
     */
    public Field createField() {
        return new Field(this);
    }


    /**
     * Gets the number of data elements contained in this field (the field-width).
     *
     * @return the field-width
     */
    public final int getNumDataElems() {
        return _numDataElems;
    }

    /**
     * Computes the size in bytes required to store the contents of a field described by this field-info and returns
     * it.
     *
     * @return the field size in bytes
     */
    public final int getSizeInBytes() {
        return getNumDataElems() * getDataTypeElemSize(getDataType());
    }

    /**
     * Returns a string representation of this field-info which can be used for debugging purposes.
     */
    @Override
    public String toString() {
        StringBuffer sb = new StringBuffer();
        sb.append("FieldInfo[");
        sb.append("'");
        sb.append(getName());
        sb.append("',");
        sb.append(ProductData.getTypeString(getDataType()));
        sb.append(",");
        sb.append(getNumDataElems());
        sb.append(",'");
        sb.append(getPhysicalUnit());
        sb.append("','");
        sb.append(getDescription());
        sb.append("']");
        return sb.toString();
    }
}
