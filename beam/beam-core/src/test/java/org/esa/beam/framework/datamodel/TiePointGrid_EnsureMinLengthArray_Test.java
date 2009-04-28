/*
 * $Id: TiePointGrid_EnsureMinLengthArray_Test.java,v 1.1 2009-04-28 14:39:33 lveci Exp $
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

public class TiePointGrid_EnsureMinLengthArray_Test extends TestCase {

    public void testEnsureMinLengthArray_ArrayIsNull() {
        final int[] ints = TiePointGrid.ensureMinLengthArray((int[]) null, 4);
        assertNotNull(ints);
        assertEquals(4, ints.length);

        final float[] floats = TiePointGrid.ensureMinLengthArray((float[]) null, 4);
        assertNotNull(floats);
        assertEquals(4, floats.length);

        final double[] doubles = TiePointGrid.ensureMinLengthArray((double[]) null, 4);
        assertNotNull(doubles);
        assertEquals(4, doubles.length);
    }

    public void testEnsureMinLengthArray_LengthIsEqual() {
        final int[] intsGiven = new int[4];
        final int[] intsReturned = TiePointGrid.ensureMinLengthArray(intsGiven, 4);
        assertSame(intsGiven, intsReturned);

        final float[] floatsGiven = new float[4];
        final float[] floatsReturned = TiePointGrid.ensureMinLengthArray(floatsGiven, 4);
        assertSame(floatsGiven, floatsReturned);

        final double[] doublesGiven = new double[4];
        final double[] doublesReturned = TiePointGrid.ensureMinLengthArray(doublesGiven, 4);
        assertSame(doublesGiven, doublesReturned);
    }

    public void testEnsureMinLengthArray_LengthIsBigger() {
        final int[] intsGiven = new int[6];
        final int[] intsReturned = TiePointGrid.ensureMinLengthArray(intsGiven, 4);
        assertSame(intsGiven, intsReturned);

        final float[] floatsGiven = new float[6];
        final float[] floatsReturned = TiePointGrid.ensureMinLengthArray(floatsGiven, 4);
        assertSame(floatsGiven, floatsReturned);

        final double[] doublesGiven = new double[6];
        final double[] doublesReturned = TiePointGrid.ensureMinLengthArray(doublesGiven, 4);
        assertSame(doublesGiven, doublesReturned);
    }

    public void testEnsureMinLengthArray_IllegalSize() {
        try {
            TiePointGrid.ensureMinLengthArray(new int[3], 4);
            fail();
        } catch (IllegalArgumentException expected) {
        } catch (Exception e) {
            fail("IllegalArgumentException expected");
        }

        try {
            TiePointGrid.ensureMinLengthArray(new float[3], 4);
            fail();
        } catch (IllegalArgumentException expected) {
        } catch (Exception e) {
            fail("IllegalArgumentException expected");
        }

        try {
            TiePointGrid.ensureMinLengthArray(new double[3], 4);
            fail();
        } catch (IllegalArgumentException expected) {
        } catch (Exception e) {
            fail("IllegalArgumentException expected");
        }
    }
}
