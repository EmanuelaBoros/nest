/*
 * $Id: ColorValidator.java,v 1.1 2009-04-28 14:39:33 lveci Exp $
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
package org.esa.beam.framework.param.validators;

import java.awt.Color;

import org.esa.beam.framework.param.AbstractParamValidator;
import org.esa.beam.framework.param.ParamConstants;
import org.esa.beam.framework.param.ParamFormatException;
import org.esa.beam.framework.param.ParamParseException;
import org.esa.beam.framework.param.ParamValidateException;
import org.esa.beam.framework.param.Parameter;
import org.esa.beam.util.Debug;
import org.esa.beam.util.StringUtils;

//@todo 1 se/** - add (more) class documentation

public class ColorValidator extends AbstractParamValidator {

    public ColorValidator() {
    }

    public Object parse(Parameter parameter, String text) throws ParamParseException {
        Debug.assertTrue(text != null);
        if (isAllowedNullText(parameter, text.trim())) {
            return null;
        }
        Color c = StringUtils.parseColor(text);
        if (c == null) {
            throw new ParamParseException(parameter, ParamConstants.ERR_MSG_INVALID_COLOR); /*I18N*/
        }
        return c;
    }

    /**
     * Returns the given {@link Color color} value as a string, according to the rules of the {@link Parameter parameter}.
     * If the value is <code>null</code> and {@link #isAllowedNullValue null value is allowed},
     * this method returns an empty string, otherwise a {@link ParamFormatException} will be thrown.
     *
     * @param parameter the parameter which contains the rules to format
     * @param value     the value to format
     *
     * @return the value as string or an empty string.
     *
     * @throws ParamFormatException if the value is <code>null</code> and
     *                              {@link #isAllowedNullValue null value is not allowed}
     *                              or the value is not an instance of {@link java.awt.Color}.
     */
    public String format(Parameter parameter, Object value) throws ParamFormatException {
        if (isAllowedNullValue(parameter, value)) {
            return "";
        }
        Color c = castToColor(value);
        if (c == null) {
            throw new ParamFormatException(parameter, ParamConstants.ERR_MSG_NOT_COLOR_TYPE);
        }
        return StringUtils.formatColor(c);
    }

    public void validate(Parameter parameter, Object value) throws ParamValidateException {
        validateThatNullValueIsAllowed(parameter, value);
        if (value == null) {
            return;
        }
        Color colorValue = castToColor(value);
        if (colorValue == null) {
            throw new ParamValidateException(parameter, ParamConstants.ERR_MSG_NOT_COLOR); /*I18N*/
        }
        validateThatValueIsInValueSet(parameter, value);
    }

    @Override
    public boolean equalValues(Parameter parameter, Object value1, Object value2) {

        Color c1 = castToColor(value1);
        Color c2 = castToColor(value2);

        if (c1 == c2) {
            return true;
        }

        if (c1 == null && c2 != null) {
            return false;
        }

        if (c1 != null && c2 == null) {
            return false;
        }

        return c1.equals(c2);
    }


    protected static Color castToColor(Object value) {
        return (value instanceof Color) ? (Color) value : null;
    }
}
