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
 * Tests for class {@link ColumnMajorMatrixFactory}.
 *
 * @author Ralf Quast
 * @version $Revision: 1.1 $ $Date: 2009-04-28 14:39:34 $
 */
public class ColumnMajorMatrixFactoryTest extends TestCase {

    private MatrixFactory matrixFactory;

    public void testCreateMatrix() {
        final double[] values = {1, 4, 2, 5, 3, 6};

        final double[][] matrix = matrixFactory.createMatrix(2, 3, values);
        assertEquals(2, matrix.length);
        assertEquals(3, matrix[0].length);
        assertEquals(3, matrix[1].length);

        assertTrue(Arrays.equals(new double[]{1, 2, 3}, matrix[0]));
        assertTrue(Arrays.equals(new double[]{4, 5, 6}, matrix[1]));
    }

    @Override
    protected void setUp() throws Exception {
        matrixFactory = new ColumnMajorMatrixFactory();
    }
}
