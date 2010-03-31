/*
 * $Id: ModisBandDescription.java,v 1.2 2010-03-31 13:59:56 lveci Exp $
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
package org.esa.beam.dataio.modis.productdb;

public class ModisBandDescription {

    private String _name;
    private boolean _isSpectral;
    private String _scaleMethod;
    private String _scaleName;
    private String _offsetName;
    private String _unitName;
    private String _bandName;
    private String _descName;
    private ModisSpectralInfo _specInfo;

    /**
     * Creates the object with given parameter set.
     *
     * @param name          the name of the band (without spectral extension)
     * @param isSpectral    whether the badnd is a spectral band or not
     * @param scalingMethod the scaling method to be used for this band (lin, exp ..)
     * @param scaleName     name of the attribute containing the scale factors
     * @param offsetName    name of the attribute containing the scale offsets
     * @param unitName      name off the attribute containing the physical unit
     * @param bandName      name of the attribute containing the spectral extensions (band names)
     * @param descName      name of the attribute containing a description of the band
     */
    public ModisBandDescription(final String name, final String isSpectral, final String scalingMethod,
                                final String scaleName, final String offsetName,
                                final String unitName, final String bandName, final String descName) {
        _name = name;
        if (isSpectral != null) {
            _isSpectral = isSpectral.equalsIgnoreCase("true");
        } else {
            // assume not
            _isSpectral = false;
        }
        _scaleMethod = scalingMethod;
        _scaleName = scaleName;
        _offsetName = offsetName;
        _unitName = unitName;
        _bandName = bandName;
        _descName = descName;
    }

    /**
     * Retrieves the name of the band
     *
     * @return the name
     */
    public String getName() {
        return _name;
    }

    public void setName(String name) {
        _name = name;
    }

    /**
     * Retrieves the scaling method to be used for this band.
     *
     * @return the scaling method
     */
    public String getScalingMethod() {
        return _scaleMethod;
    }

    /**
     * Retrieves the name of the attribute containing the scaling offset
     *
     * @return the name of the attribute
     */
    public String getOffsetAttribName() {
        return _offsetName;
    }

    /**
     * Retrieves the name of the attribute containing the scaling factor
     *
     * @return the name of the attribute
     */
    public String getScaleAttribName() {
        return _scaleName;
    }

    /**
     * Retrieves the name of the attribute containing the physical unit
     *
     * @return the name of the attribute
     */
    public String getUnitAttribName() {
        return _unitName;
    }

    /**
     * Retrieves the name of the attribute containing the band names (spectral extensions)
     *
     * @return the name of the attribute
     */
    public String getBandAttribName() {
        return _bandName;
    }

    /**
     * Retrieves whether this band is a spectral band or not
     *
     * @return <code>true</code> if this band is a spectral band, otherwise <code>false</code>.
     */
    public boolean isSpectral() {
        return _isSpectral;
    }

    /**
     * Retrieves the name of the attribute containing the band descritpion.
     *
     * @return the name of the attribute
     */
    public String getDescriptionAttribName() {
        return _descName;
    }

    public void setSpecInfo(final ModisSpectralInfo specInfo) {
        _specInfo = specInfo;
    }

    public ModisSpectralInfo getSpecInfo() {
        return _specInfo;
    }

    public boolean hasSpectralInfo() {
        return _specInfo != null;
    }
}
