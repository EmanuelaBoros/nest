/*
 * $Id: UnconstrainedLSU.java,v 1.1 2009-04-28 14:39:33 lveci Exp $
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
package org.esa.beam.util.math;

import Jama.Matrix;

import java.text.MessageFormat;

/**
 * Performs an unconstrained linear spectral unmixing.
 *
 * @author Ralf Quast
 * @author Helmut Schiller (GKSS)
 * @version $Revision: 1.1 $ $Date: 2009-04-28 14:39:33 $
 * @since 4.1
 */
public class UnconstrainedLSU implements SpectralUnmixing {

    private final double[][] endmemberMatrix;
    private final double[][] inverseEndmemberMatrix;

    /**
     * Constructs a new instance of this class.
     *
     * @param endmembers the endmembers, where
     *                   number of rows = number of spectral channels
     *                   number of cols = number of endmember spectra
     */
    public UnconstrainedLSU(double[][] endmembers) {
        if (!LinearAlgebra.isMatrix(endmembers)) {
            throw new IllegalArgumentException("Parameter 'endmembers' is not a matrix.");
        }

        endmemberMatrix = endmembers;
        inverseEndmemberMatrix = new Matrix(endmembers).inverse().getArrayCopy();
    }

    /**
     * Returns the endmembers.
     *
     * @return endmembers the endmembers, where
     *         number of rows = number of spectral channels
     *         number of cols = number of endmember spectra
     */
    public double[][] getEndmembers() {
        return endmemberMatrix;
    }

    @Override
    public double[][] unmix(double[][] spectra) {
        final int actualRowCount = spectra.length;
        final int expectedRowCount = endmemberMatrix.length;

        if (actualRowCount != expectedRowCount) {
            throw new IllegalArgumentException(MessageFormat.format(
                    "Parameter ''spectra'' is not a matrix with {0} rows.", expectedRowCount));
        }

        return LinearAlgebra.multiply(inverseEndmemberMatrix, spectra);
    }

    @Override
    public double[][] mix(double[][] abundances) {
        final int actualRowCount = abundances.length;
        final int expectedRowCount = endmemberMatrix[0].length;

        if (actualRowCount != expectedRowCount) {
            throw new IllegalArgumentException(MessageFormat.format(
                    "Parameter ''abundances'' is not a matrix with {0} rows.", expectedRowCount));
        }

        return LinearAlgebra.multiply(endmemberMatrix, abundances);
    }
}
