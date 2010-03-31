/*
 * $Id: ApproximatorTest.java,v 1.2 2010-03-31 13:59:56 lveci Exp $
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

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

public class ApproximatorTest extends TestCase {

    public final static double EPS = 1e-10;
    public final static int N = 72 * 72;

    public ApproximatorTest(String s) {
        super(s);
    }

    public static Test suite() {
        return new TestSuite(ApproximatorTest.class);
    }

    public void testFX() {
        final double axx = -3.0;
        final double ax = 2.0;
        final double a1 = -1.0;

        double data[][] = new double[N][2];
        for (int i = 0; i < data.length; i++) {
            double x = random(-8, 8);
            double y = axx * x * x + ax * x + a1;
            data[i][0] = x;
            data[i][1] = y;
        }
        final FX[] funcs = new FX[]{FX.XX, FX.X, FX.ONE};
        double[] coeffs = new double[3];
        Approximator.approximateFX(data, null, funcs, coeffs);
        double rmse = Approximator.getRMSE(data, null, funcs, coeffs);

        assertEquals(axx, coeffs[0], EPS);
        assertEquals(ax, coeffs[1], EPS);
        assertEquals(a1, coeffs[2], EPS);
        assertEquals(0.0, rmse, EPS);
    }

    public void testFXY() {
        final double axx = -6.0;
        final double ayy = +5.0;
        final double axy = -4.0;
        final double ax = +3.0;
        final double ay = -2.0;
        final double a1 = +1.0;

        double data[][] = new double[N][3];
        for (int i = 0; i < data.length; i++) {
            double x = random(-8, 8);
            double y = random(-8, 8);
            double z = axx * x * x + ayy * y * y + axy * x * y + ax * x + ay * y + a1;
            data[i][0] = x;
            data[i][1] = y;
            data[i][2] = z;
        }
        final FXY[] funcs = new FXY[]{FXY.X2, FXY.Y2, FXY.XY,
                                      FXY.X, FXY.Y,
                                      FXY.ONE};
        double[] coeffs = new double[6];
        Approximator.approximateFXY(data, null, funcs, coeffs);
        double rmse = Approximator.computeRMSE(data, null, funcs, coeffs);

        assertEquals(axx, coeffs[0], EPS);
        assertEquals(ayy, coeffs[1], EPS);
        assertEquals(axy, coeffs[2], EPS);
        assertEquals(ax, coeffs[3], EPS);
        assertEquals(ay, coeffs[4], EPS);
        assertEquals(a1, coeffs[5], EPS);
        assertEquals(0.0, rmse, EPS);
    }

    private static double random(double x1, double x2) {
        return x1 + Math.random() * (x2 - x1);
    }
}
