/*
 * $Id: ModisTiePointDescription.java,v 1.1 2009-09-25 19:03:49 lveci Exp $
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
package org.esa.beam.dataio.modis.productdb;

public class ModisTiePointDescription {

    private String _name;
    private String _scale;
    private String _offset;
    private String _unit;

    /**
     * Creates the object with given parameters.
     *
     * @param name            the name of the tie point grid
     * @param scaleAttribute  the name of the attribute containing the scale factor
     * @param offsetAttribute the name of the attribute containing the scaling offset
     * @param unitAttribute   the name of the attribute containi8ng the unit name
     */
    public ModisTiePointDescription(final String name, final String scaleAttribute, final String offsetAttribute,
                                    final String unitAttribute) {
        _name = name;
        _scale = scaleAttribute;
        _offset = offsetAttribute;
        _unit = unitAttribute;
    }

    /**
     * Retrievs the name of the rtie point grid
     *
     * @return the name
     */
    public String getName() {
        return _name;
    }

    /**
     * Retrieves the name of the scaling factor attribute
     *
     * @return the name
     */
    public String getScaleAttribName() {
        return _scale;
    }

    /**
     * Retrieves the name of the scaling offset attribute
     *
     * @return the name
     */
    public String getOffsetAttribName() {
        return _offset;
    }

    /**
     * Retrieves the name of the attribute for the physical unit
     *
     * @return the name
     */
    public String getUnitAttribName() {
        return _unit;
    }
}
