/*
 * $Id: LineFigure.java,v 1.2 2010-02-08 21:57:50 lveci Exp $
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
package org.esa.beam.framework.draw;

import java.awt.Shape;
import java.util.Map;

/**
 * A shape figure that represents a one-dimensional line.
 *
 * @deprecated since BEAM 4.7, no replacement
 */
@Deprecated
public class LineFigure extends ShapeFigure {

    public LineFigure(Shape shape, Map<String, Object> attributes) {
        super(shape, true, attributes);
    }

}
