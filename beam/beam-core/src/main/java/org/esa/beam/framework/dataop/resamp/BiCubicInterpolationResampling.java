/*
 * Copyright (C) 2010 Brockmann Consult GmbH (info@brockmann-consult.de)
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation; either version 3 of the License, or (at your option)
 * any later version.
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, see http://www.gnu.org/licenses/
 */

package org.esa.beam.framework.dataop.resamp;


import javax.xml.crypto.dsig.keyinfo.KeyInfo;

final class BiCubicInterpolationResampling implements Resampling {

    private final static float[][] invA = {
                {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                {0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                {-3, 3, 0, 0, -2, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                {2, -2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                {0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0},
                {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0},
                {0, 0, 0, 0, 0, 0, 0, 0, -3, 3, 0, 0, -2, -1, 0, 0},
                {0, 0, 0, 0, 0, 0, 0, 0, 2, -2, 0, 0, 1, 1, 0, 0},
                {-3, 0, 3, 0, 0, 0, 0, 0, -2, 0, -1, 0, 0, 0, 0, 0},
                {0, 0, 0, 0, -3, 0, 3, 0, 0, 0, 0, 0, -2, 0, -1, 0},
                {9, -9, -9, 9, 6, 3, -6, -3, 6, -6, 3, -3, 4, 2, 2, 1},
                {-6, 6, 6, -6, -3, -3, 3, 3, -4, 4, -2, 2, -2, -2, -1, -1},
                {2, 0, -2, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0},
                {0, 0, 0, 0, 2, 0, -2, 0, 0, 0, 0, 0, 1, 0, 1, 0},
                {-6, 6, 6, -6, -4, -2, 4, 2, -3, 3, -3, 3, -2, -1, -2, -1},
                {4, -4, -4, 4, 2, 2, -2, -2, 2, -2, 2, -2, 1, 1, 1, 1}};

    public String getName() {
        return "BICUBIC_INTERPOLATION";
    }

    public final Index createIndex() {
        return new Index(2, 1);
    }

    public final void computeIndex(final float x,
                                   final float y,
                                   final int width,
                                   final int height,
                                   final Index index) {
        index.x = x;
        index.y = y;
        index.width = width;
        index.height = height;

        final int i0 = (int) Math.floor(x);
        final int j0 = (int) Math.floor(y);


        final float di = x - (i0 + 0.5f);
        final float dj = y - (j0 + 0.5f);

        index.i0 = i0;
        index.j0 = j0;

        final int iMax = width - 1;
        final int jMax = height - 1;

        if (di >= 0) {
            final int i1 = i0 + 1;
            index.i[0] = (i0 < 0) ? 0 : (i0 > iMax) ? iMax : i0; //Index.crop(i0, iMax);
            index.i[1] = (i1 < 0) ? 0 : (i1 > iMax) ? iMax : i1; //Index.crop(i0 + 1, iMax);
            index.ki[0] = di;
        } else {
            final int i1 = i0 - 1;
            index.i[0] = (i1 < 0) ? 0 : (i1 > iMax) ? iMax : i1; //Index.crop(i0 - 1, iMax);
            index.i[1] = (i0 < 0) ? 0 : (i0 > iMax) ? iMax : i0; //Index.crop(i0, iMax);
            index.ki[0] = di + 1;
        }

        if (dj >= 0) {
            final int j1 = j0 + 1;
            index.j[0] = (j0 < 0) ? 0 : (j0 > jMax) ? jMax : j0; //Index.crop(j0, jMax);
            index.j[1] = (j1 < 0) ? 0 : (j1 > jMax) ? jMax : j1; //Index.crop(j0 + 1, jMax);
            index.kj[0] = dj;
        } else {
            final int j1 = j0 - 1;
            index.j[0] = (j1 < 0) ? 0 : (j1 > jMax) ? jMax : j1; //Index.crop(j0 - 1, jMax);
            index.j[1] = (j0 < 0) ? 0 : (j0 > jMax) ? jMax : j0; //Index.crop(j0, jMax);
            index.kj[0] = dj + 1;
        }

        /*
        final float di = x - i0;
        final float dj = y - j0;

        index.i0 = i0;
        index.j0 = j0;

        final int iMax = width - 1;
        final int jMax = height - 1;

        index.i[0] = Index.crop(i0, iMax);
        index.i[1] = Index.crop(i0 + 1, iMax);
        index.ki[0] = di;

        index.j[0] = Index.crop(j0, jMax);
        index.j[1] = Index.crop(j0 + 1, jMax);
        index.kj[0] = dj;
        */
    }

    public final float resample(final Raster raster,
                                final Index index) throws Exception {

        final float[][] v = new float[4][4];
        final int x0 = index.i[0] - 1;
        final int y0 = index.j[0] - 1;
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                final int x = Index.crop(x0 + i, index.width-1);
                final int y = Index.crop(y0 + j, index.height-1);
                v[j][i] = raster.getSample(x, y);
                if(Double.isNaN(v[j][i])) {
                    return raster.getSample(index.i0, index.j0);
                }
            }
        }

        // the four grid points of a rectangular grid cell are numbered as the following:
        // p1    p2
        //
        // p3    p4

        final float[] z = new float[4];     // function values
        final float[] z1 = new float[4];    // 1st order derivative in y direction
        final float[] z2 = new float[4];    // 1st order derivative in x direction
        final float[] z12 = new float[4];   // cross derivative

        z[0] = v[1][1];
        z[1] = v[1][2];
        z[2] = v[2][1];
        z[3] = v[2][2];

        z1[0] = (v[1][2] - v[1][0]) / 2.0f;
        z1[1] = (v[1][3] - v[1][1]) / 2.0f;
        z1[2] = (v[2][2] - v[2][0]) / 2.0f;
        z1[3] = (v[2][3] - v[2][1]) / 2.0f;

        z2[0] = (v[2][1] - v[0][1]) / 2.0f;
        z2[1] = (v[2][2] - v[0][2]) / 2.0f;
        z2[2] = (v[3][1] - v[1][1]) / 2.0f;
        z2[3] = (v[3][2] - v[1][2]) / 2.0f;

        z12[0] = (v[2][2] - v[2][0] - v[0][2] + v[0][0]) / 4.0f;
        z12[1] = (v[2][3] - v[2][1] - v[0][3] + v[0][1]) / 4.0f;
        z12[2] = (v[3][2] - v[3][0] - v[1][2] + v[1][0]) / 4.0f;
        z12[3] = (v[3][3] - v[3][1] - v[1][3] + v[1][1]) / 4.0f;

        return bcuint(z, z1, z2, z12, index.ki[0], index.kj[0]);
    }

	private static float bcuint(final float z[], final float z1[], final float z2[],
                         final float z12[], final float t, final float u) {

        // alpha = [a00 a10 a20 a30 a01 a11 a21 a31 a02 a12 a22 a32 a03 a13 a23 a33]
		final float[][] a = new float[4][4];
		bcucof(z, z1, z2, z12, a);

		float ansy = 0.0f;
		for (int i = 3; i >= 0; i--) {
			ansy = t*ansy + ((a[i][3]*u + a[i][2])*u + a[i][1])*u + a[i][0];
		}

        // todo v is not used?? check with Jun
        /*
        float t2 = t*t;
        float t3 = t2*t;
        float u2 = u*u;
        float u3 = u2*u;

        float v = a[0][0] + a[0][1]*u + a[0][2]*u2 + a[0][3]*u3 +
                  a[1][0]*t + a[1][1]*t*u + a[1][2]*t*u2 + a[1][3]*t*u3 +
                  a[2][0]*t2 + a[2][1]*t2*u + a[2][2]*t2*u2 + a[2][3]*t2*u3 +
                  a[3][0]*t3 + a[3][1]*t3*u + a[3][2]*t3*u2 + a[3][3]*t3*u3;
        */
		return ansy;
	}

	private static void bcucof(final float z[], final float z1[], final float z2[], final float z12[],
                               final float[][] a) {

        // x = [f(0,0) f(1,0) f(0,1) f(1,1) fx(0,0) fx(1,0) fx(0,1) fx(1,1) fy(0,0) fy(1,0) fy(0,1) fy(1,1) fxy(0,0) fxy(1,0) fxy(0,1) fxy(1,1)]
        // alpha = [a00 a10 a20 a30 a01 a11 a21 a31 a02 a12 a22 a32 a03 a13 a23 a33]
        // alpha = invA*x

        final float[] x = new float[16];
        for (int i = 0; i < 4; i++) {
            x[i] = z[i];
            x[i+4] = z1[i];
            x[i+8] = z2[i];
            x[i+12] = z12[i];
        }

		final float[] cl = new float[16];
		for (int i = 0; i < 16; i++) {
			float xx = 0.0f;
			for (int k = 0; k < 16; k++) {
                xx += invA[i][k]*x[k];
			}
			cl[i] = xx;
		}
		
		int l = 0;
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				a[j][i] = cl[l++];
			}
		}
	}

    @Override
    public String toString() {
        return "BiCubic interpolation resampling";
    }
}
