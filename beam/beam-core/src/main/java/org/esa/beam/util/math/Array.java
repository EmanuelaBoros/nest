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

/**
 * Interface for wrapping primitive arrays.
 *
 * @author Ralf Quast
 * @version $Revision: 1.1 $ $Date: 2009-04-28 14:39:33 $
 */
public interface Array {

    /**
     * Returns the length of the primitive array wrapped.
     *
     * @return the lenght of the primitive array wrapped.
     */
    public abstract int getLength();

    /**
     * Returns the ith value of the array wrapped.
     *
     * @param i the array index.
     * @return the ith value of the array wrapped.
     */
    public abstract double getValue(int i);

    /**
     * Copies an array from the primitive array wrapped, beginning at the
     * specified position, to the specified position of the destination array.
     *
     * @param srcPos  starting position in the primitive array wrapped.
     * @param dest    the destination array.
     * @param destPos starting position in the destination array.
     * @param length  the number of array elements to be copied.
     */
    public abstract void copyTo(int srcPos, double[] dest, int destPos, int length);

    /**
     * Class for wrapping {@code double} primitive arrays.
     */
    public static class Double implements Array {
        private final double[] values;

        public Double(double... values) {
            if (values == null) {
                throw new NullPointerException("values == null");
            }
            this.values = values;
        }

        @Override
        public final int getLength() {
            return values.length;
        }

        @Override
        public final double getValue(int i) {
            return values[i];
        }

        @Override
        public final void copyTo(int srcPos, double[] dest, int destPos, int length) {
            System.arraycopy(values, srcPos, dest, destPos, length);
        }
    }

    /**
     * Class for wrapping {@code float} primitive arrays.
     */
    public static class Float implements Array {
        private final float[] values;

        public Float(float... values) {
            if (values == null) {
                throw new IllegalArgumentException("values == null");
            }
            this.values = values;
        }

        @Override
        public final int getLength() {
            return values.length;
        }

        @Override
        public final double getValue(int i) {
            return values[i];
        }

        @Override
        public final void copyTo(int srcPos, double[] dest, int destPos, int length) {
            for (int i = 0; i < length; ++i) {
                dest[destPos + i] = values[srcPos + i];
            }
        }
    }
}
