/* 
 * Copyright (C) 2002-2008 by Brockmann Consult
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
package org.esa.beam.util.math;

import junit.framework.TestCase;

import java.util.Arrays;

/**
 * Tests for class {@link VectorLookupTable}.
 *
 * @author Ralf Quast
 * @version $Revision: 1.2 $ $Date: 2010-03-31 13:59:56 $
 */
public class VectorLookupTableTest extends TestCase {

    public void testArrayInterpolation1D() {
        final double[] dimension = new double[]{0, 1};
        final double[] values = new double[]{0, 2, 1, 3};

        final VectorLookupTable lut = new VectorLookupTable(2, values, dimension);
        assertEquals(1, lut.getDimensionCount());

        assertEquals(0.0, lut.getDimension(0).getMin(), 0.0);
        assertEquals(1.0, lut.getDimension(0).getMax(), 0.0);

        assertTrue(Arrays.equals(new double[]{0.0, 2.0}, lut.getValues(0.0)));
        assertTrue(Arrays.equals(new double[]{1.0, 3.0}, lut.getValues(1.0)));
        assertTrue(Arrays.equals(new double[]{0.5, 2.5}, lut.getValues(0.5)));
    }

    public void testArrayInterpolation2D() {
        final double[][] dimensions = new double[][]{{0, 1}, {0, 1}};
        final double[] values = new double[]{0, 4, 1, 5, 2, 6, 3, 7};

        final VectorLookupTable lut = new VectorLookupTable(2, values, dimensions);
        assertEquals(2, lut.getDimensionCount());

        assertEquals(0.0, lut.getDimension(0).getMin(), 0.0);
        assertEquals(1.0, lut.getDimension(0).getMax(), 0.0);
        assertEquals(0.0, lut.getDimension(1).getMin(), 0.0);
        assertEquals(1.0, lut.getDimension(1).getMax(), 0.0);

        assertTrue(Arrays.equals(new double[]{0.0, 4.0}, lut.getValues(0.0, 0.0)));
        assertTrue(Arrays.equals(new double[]{1.0, 5.0}, lut.getValues(0.0, 1.0)));
        assertTrue(Arrays.equals(new double[]{2.0, 6.0}, lut.getValues(1.0, 0.0)));
        assertTrue(Arrays.equals(new double[]{3.0, 7.0}, lut.getValues(1.0, 1.0)));

        assertTrue(Arrays.equals(new double[]{0.5, 4.5}, lut.getValues(0.0, 0.5)));
        assertTrue(Arrays.equals(new double[]{1.5, 5.5}, lut.getValues(0.5, 0.5)));
        assertTrue(Arrays.equals(new double[]{2.5, 6.5}, lut.getValues(1.0, 0.5)));
    }

    public void testArrayInterpolationWithFloatValues() {
        final float[][] dimensions = new float[][]{{0, 1}, {0, 1}};
        final float[] values = new float[]{0, 4, 1, 5, 2, 6, 3, 7};

        final VectorLookupTable lut = new VectorLookupTable(2, values, dimensions);
        assertTrue(Arrays.equals(new double[]{0.0, 4.0}, lut.getValues(0.0, 0.0)));
        assertTrue(Arrays.equals(new double[]{1.0, 5.0}, lut.getValues(0.0, 1.0)));
        assertTrue(Arrays.equals(new double[]{2.0, 6.0}, lut.getValues(1.0, 0.0)));
        assertTrue(Arrays.equals(new double[]{3.0, 7.0}, lut.getValues(1.0, 1.0)));

        assertTrue(Arrays.equals(new double[]{0.5, 4.5}, lut.getValues(0.0, 0.5)));
        assertTrue(Arrays.equals(new double[]{1.5, 5.5}, lut.getValues(0.5, 0.5)));
        assertTrue(Arrays.equals(new double[]{2.5, 6.5}, lut.getValues(1.0, 0.5)));
    }
}
