/*
 * $Id: HdfUtils.java,v 1.1 2009-09-25 19:03:49 lveci Exp $
 *
 * Copyright (C) 2002,2003  by Brockmann Consult (info@brockmann-consult.de)
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation. This program is distributed in the hope it will
 * be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 */
package org.esa.beam.dataio.modis.hdf;

import ncsa.hdf.hdflib.HDFConstants;
import ncsa.hdf.hdflib.HDFException;
import ncsa.hdf.hdflib.HDFNativeData;
import org.esa.beam.dataio.modis.hdf.lib.HDF;
import org.esa.beam.framework.datamodel.MetadataAttribute;
import org.esa.beam.framework.datamodel.ProductData;
import org.esa.beam.util.Debug;
import org.esa.beam.util.logging.BeamLogManager;

import java.util.ArrayList;
import java.util.List;

public class HdfUtils {

    /**
     * Decodes the hdf data type into a product data type.
     *
     * @param hdfType
     *
     * @return product data type
     *
     * @see ProductData
     */
    public static int decodeHdfDataType(final int hdfType) {
        switch (hdfType) {
        case HDFConstants.DFNT_UCHAR8:
        case HDFConstants.DFNT_UINT8:
            return ProductData.TYPE_UINT8;

        case HDFConstants.DFNT_CHAR8:
        case HDFConstants.DFNT_INT8:
            return ProductData.TYPE_INT8;

        case HDFConstants.DFNT_INT16:
            return ProductData.TYPE_INT16;

        case HDFConstants.DFNT_UINT16:
            return ProductData.TYPE_UINT16;

        case HDFConstants.DFNT_INT32:
            return ProductData.TYPE_INT32;

        case HDFConstants.DFNT_UINT32:
            return ProductData.TYPE_UINT32;

        case HDFConstants.DFNT_FLOAT32:
            return ProductData.TYPE_FLOAT32;

        case HDFConstants.DFNT_FLOAT64:
            return ProductData.TYPE_FLOAT64;

        default:
            return ProductData.TYPE_UNDEFINED;
        }
    }

    /**
     * Decodes a byte buffer as read from the HDF library to a useable string.
     *
     * @return a new {@link HdfAttributeContainer} with the decoded string as value.
     */
    public static HdfAttributeContainer decodeByteBufferToAttribute(byte[] buf, int NT, int count, String name)
            throws HDFException {
        final String strVal;
        final int incr = HDF.getWrap().DFKNTsize(NT);

        String strTemp;

        switch (NT) {
        case HDFConstants.DFNT_CHAR:
        case HDFConstants.DFNT_UCHAR8:
            strVal = new String(buf, 0, count).trim();
            break;
        case HDFConstants.DFNT_UINT8:
            strTemp = "";
            for (int n = 0; n < count; n++) {
                Byte bval = buf[n];
                Short shval;
                if (bval.shortValue() < 0) {
                    shval = (short) (bval.intValue() + 256);
                } else {
                    shval = bval.shortValue();
                }
                strTemp += shval.toString();
                strTemp += (",");
            }
            strVal = strTemp.substring(0, strTemp.length() - 1);
            break;
        case HDFConstants.DFNT_INT8:
            strTemp = "";
            for (int n = 0; n < count; n++) {
                final int pos = n * incr;
                Byte bval = buf[pos];
                strTemp += bval.toString();
                strTemp += (",");
            }
            strVal = strTemp.substring(0, strTemp.length() - 1);
            break;
        case HDFConstants.DFNT_INT16:
            strTemp = "";
            for (int n = 0; n < count; n++) {
                final int pos = n * incr;
                Short shval = HDFNativeData.byteToShort(buf, pos);
                strTemp += shval.toString();
                strTemp += (",");
            }
            strVal = strTemp.substring(0, strTemp.length() - 1);
            break;
        case HDFConstants.DFNT_UINT16:
            strTemp = "";
            for (int n = 0; n < count; n++) {
                final int pos = n * incr;
                Short shval = HDFNativeData.byteToShort(buf, pos);
                Integer ival;
                if (shval < 0) {
                    ival = (shval.intValue() + 65536);
                } else {
                    ival = shval.intValue();
                }
                strTemp += ival.toString();
                strTemp += (",");
            }
            strVal = strTemp.substring(0, strTemp.length() - 1);
            break;
        case HDFConstants.DFNT_INT32:
            strTemp = "";
            for (int n = 0; n < count; n++) {
                final int pos = n * incr;
                Integer ival = HDFNativeData.byteToInt(buf, pos);
                strTemp += ival.toString();
                strTemp += (",");
            }
            strVal = strTemp.substring(0, strTemp.length() - 1);
            break;
        case HDFConstants.DFNT_UINT32:
            strTemp = "";
            for (int n = 0; n < count; n += incr) {
                final int pos = n * incr;
                Integer ival = HDFNativeData.byteToInt(buf, pos);
                Long lVal;

                if (ival < 0) {
                    lVal = ival.longValue() + 4294967296L;
                } else {
                    lVal = ival.longValue();
                }
                strTemp += lVal.toString();
                strTemp += (",");
            }
            strVal = strTemp.substring(0, strTemp.length() - 1);
            break;
        case HDFConstants.DFNT_FLOAT32:
            strTemp = "";
            for (int n = 0; n < count; n++) {
                final int pos = n * incr;
                Float fval = HDFNativeData.byteToFloat(buf, pos);
                strTemp += fval.toString();
                strTemp += (",");
            }
            strVal = strTemp.substring(0, strTemp.length() - 1);
            break;
        case HDFConstants.DFNT_DOUBLE:
            strTemp = "";
            for (int n = 0; n < count; n++) {
                final int pos = n * incr;
                Double dval = HDFNativeData.byteToDouble(buf, pos);
                strTemp += dval.toString();
                strTemp += (",");
            }
            strVal = strTemp.substring(0, strTemp.length() - 1);
            break;
        default:
            return null;
        }
        return new HdfAttributeContainer(name, NT, strVal, count);
    }

    /**
     * Converts a hdf attribute container to a product metadata attribute
     *
     * @param container
     *
     * @return the product metadata attribute
     */
    public static MetadataAttribute attributeToMetadata(HdfAttributeContainer container) {
        ProductData prodData = null;
        MetadataAttribute attrib = null;

        switch (container.getHdfType()) {
        case HDFConstants.DFNT_CHAR:
        case HDFConstants.DFNT_UCHAR8:
            prodData = ProductData.createInstance(container.getStringValue());
            break;

        case HDFConstants.DFNT_UINT8:
        case HDFConstants.DFNT_INT8:
        case HDFConstants.DFNT_UINT16:
        case HDFConstants.DFNT_INT16:
        case HDFConstants.DFNT_INT32:
        case HDFConstants.DFNT_UINT32:
            prodData = ProductData.createInstance(container.getIntValues());
            break;

        case HDFConstants.DFNT_FLOAT32:
            prodData = ProductData.createInstance(container.getFloatValues());
            break;

        case HDFConstants.DFNT_DOUBLE:
            prodData = ProductData.createInstance(container.getDoubleValues());
            break;
        }

        if (prodData != null) {
            attrib = new MetadataAttribute(container.getName(), prodData, true);
        }
        return attrib;
    }

    /**
     * Reads the global attributes from the hdf file passed in.
     *
     * @param sdId the HD interface identifier of the file
     *
     * @return an instance of HdfGlobalAttributes which contains all the
     *         global hdf attributes read from file.
     *
     * @throws ncsa.hdf.hdflib.HDFException -
     */
    final static public HdfAttributes readAttributes(final int sdId) throws HDFException {
        Debug.trace("reading global attributes ...");

        final List<HdfAttributeContainer> attributes = new ArrayList<HdfAttributeContainer>();

        final int[] fileInfo = new int[2];
        // request number of datasets (fileInfo[0]) and number of global attributes (fileInfo[1])
        if (HDF.getWrap().SDfileinfo(sdId, fileInfo)) {
            final int[] sdAttrInfo = new int[2];
            final String[] sdVal = new String[1];
            for (int n = 0; n < fileInfo[1]; n++) {
                sdVal[0] = "";
                if (HDF.getWrap().SDattrinfo(sdId, n, sdVal, sdAttrInfo)) {
                    final int attrSize = HDF.getWrap().DFKNTsize(sdAttrInfo[0]) * sdAttrInfo[1] + 1;
                    final byte[] buf = new byte[attrSize];

                    if (HDF.getWrap().SDreadattr(sdId, n, buf)) {
                        final String attrName = sdVal[0].trim();
                        final HdfAttributeContainer attribute
                                = decodeByteBufferToAttribute(buf, sdAttrInfo[0], sdAttrInfo[1], attrName);

                        if (attribute != null) {
                            attributes.add(attribute);
                            Debug.trace("... " + attrName + ": " + attribute.getStringValue());
                        }
                    }
                }
            }
        }

        Debug.trace("... success");
        return new HdfAttributes(attributes);
    }

    /**
     * Retrieves a string attribute with given name from the sds specified.
     *
     * @@param sdsId
     * @@param name
     * @@return a string attribute
     */
    public static String getNamedStringAttribute(int sdsId, String name) throws HDFException {
        if (name == null) {
            return null;
        }

        final int attrIdx = HDF.getWrap().SDfindattr(sdsId, name);
        if (attrIdx == HDFConstants.FAIL) {
            return null;
        }

        final int[] attrInfo = new int[2];
        final String[] dsName = new String[]{""};
        if (HDF.getWrap().SDattrinfo(sdsId, attrIdx, dsName, attrInfo)) {
            final int attrSize = HDF.getWrap().DFKNTsize(attrInfo[0]) * attrInfo[1];
            final byte[] buf = new byte[attrSize];
            if (HDF.getWrap().SDreadattr(sdsId, attrIdx, buf)) {
                return new String(buf).trim();
            }
        }
        BeamLogManager.getSystemLogger().warning("Unable to access the attribute '" + name + '\'');
        return null;
    }
}
