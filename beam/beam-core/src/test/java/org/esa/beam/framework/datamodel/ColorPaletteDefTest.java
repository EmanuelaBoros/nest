/*
 * $Id: ColorPaletteDefTest.java,v 1.2 2010-03-31 13:59:56 lveci Exp $
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

import junit.framework.TestCase;

import java.awt.Color;

public class ColorPaletteDefTest extends TestCase {

    public void testConstructors() {
        ColorPaletteDef cpd = new ColorPaletteDef(-1.0, 1.0);
        assertEquals(256, cpd.getNumColors());
        assertEquals(3, cpd.getNumPoints());
        assertEquals(-1.0, cpd.getPointAt(0).getSample(), 1e-10);
        assertEquals(+0.0, cpd.getPointAt(1).getSample(), 1e-10);
        assertEquals(+1.0, cpd.getPointAt(2).getSample(), 1e-10);
        assertEquals(Color.BLACK, cpd.getPointAt(0).getColor());
        assertEquals(Color.GRAY, cpd.getPointAt(1).getColor());
        assertEquals(Color.WHITE, cpd.getPointAt(2).getColor());

        cpd = new ColorPaletteDef(-1.0, 0.5, 1.0);
        assertEquals(256, cpd.getNumColors());
        assertEquals(3, cpd.getNumPoints());
        assertEquals(-1.0, cpd.getPointAt(0).getSample(), 1e-10);
        assertEquals(+0.5, cpd.getPointAt(1).getSample(), 1e-10);
        assertEquals(+1.0, cpd.getPointAt(2).getSample(), 1e-10);
        assertEquals(Color.BLACK, cpd.getPointAt(0).getColor());
        assertEquals(Color.GRAY, cpd.getPointAt(1).getColor());
        assertEquals(Color.WHITE, cpd.getPointAt(2).getColor());

        cpd = new ColorPaletteDef(new ColorPaletteDef.Point[]{
                new ColorPaletteDef.Point(100, Color.ORANGE),
                new ColorPaletteDef.Point(200, Color.MAGENTA),
                new ColorPaletteDef.Point(500, Color.BLUE),
                new ColorPaletteDef.Point(600, Color.WHITE)
        });
        assertEquals(4, cpd.getNumPoints());
        assertEquals(256, cpd.getNumColors());


        cpd = new ColorPaletteDef(new ColorPaletteDef.Point[]{
                new ColorPaletteDef.Point(100, Color.ORANGE),
                new ColorPaletteDef.Point(200, Color.MAGENTA),
                new ColorPaletteDef.Point(500, Color.BLUE),
                new ColorPaletteDef.Point(600, Color.WHITE)
        }, 512);
        assertEquals(4, cpd.getNumPoints());
        assertEquals(512, cpd.getNumColors());
    }



    public void testPaletteCreation() {

        ColorPaletteDef cpd = new ColorPaletteDef(new ColorPaletteDef.Point[]{
                new ColorPaletteDef.Point(100, Color.ORANGE),
                new ColorPaletteDef.Point(200, Color.MAGENTA),
                new ColorPaletteDef.Point(300, Color.RED),
                new ColorPaletteDef.Point(400, Color.GREEN),
        }, 4);
        Color[] palette = cpd.createColorPalette(Scaling.IDENTITY);
        assertNotNull(palette);
        assertEquals(4, palette.length);
        assertEquals(Color.ORANGE, palette[0]);
        assertEquals(Color.MAGENTA, palette[1]);
        assertEquals(Color.RED, palette[2]);
        assertEquals(Color.GREEN, palette[3]);

        cpd = new ColorPaletteDef(new ColorPaletteDef.Point[]{
                new ColorPaletteDef.Point(100, Color.WHITE),
                new ColorPaletteDef.Point(200, Color.BLUE),
                new ColorPaletteDef.Point(300, Color.RED),
                new ColorPaletteDef.Point(400, Color.GREEN),
        }, 7);
        palette = cpd.createColorPalette(Scaling.IDENTITY);
        assertNotNull(palette);
        assertEquals(7, palette.length);
        assertEquals(new Color(255, 255, 255), palette[0]);
        assertEquals(new Color(128, 128, 255), palette[1]);
        assertEquals(new Color(0, 0, 255), palette[2]);
        assertEquals(new Color(128, 0, 128), palette[3]);
        assertEquals(new Color(255, 0, 0), palette[4]);
        assertEquals(new Color(128, 128, 0), palette[5]);
        assertEquals(new Color(0, 255, 0), palette[6]);
    }
}