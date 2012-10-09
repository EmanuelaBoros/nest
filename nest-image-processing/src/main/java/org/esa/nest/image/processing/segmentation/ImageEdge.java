package org.esa.nest.image.processing.segmentation;

import java.util.ArrayList;
import ij.*;
import ij.process.*;
import ij.plugin.filter.RankFilters;
import no.uib.cipr.matrix.Matrices;
import no.uib.cipr.matrix.Matrix;

/**
 * Class containing some static methods for Edge Detection, using median filter,
 * Canny-Deriche filtering, double tresholding, hysteresis tresholding. These
 * are meant to be implemented in another plugin, hence the basic functionality.
 * Methods based on the plugins deriche_ and hysteresis_
 *
 * @author thomas.boudier@snv.jussieu.fr
 * @author Joris FA Meys @created	4/2/2009
 * @version 1.0
 *
 */
public class ImageEdge {

    static final int BYTE = 0, SHORT = 1, FLOAT = 2, OTHER = 3;
    static private boolean typed = false;
    static private int type;

    /**
     * Edge detection by subsequent application of a median filter, a
     * Canny-Deriche filter double tresholding and hysteresis. This method works
     * especially well for detection of area edges in pictures with a high
     * granularity, like confocal pics at high magnification. Good results were
     * obtained with radius = 5, alpha = 0.5, upper = 100, lower = 50
     *
     * @param ip	ImageProcessor
     * @param radius	radius for the median filter
     * @param alpha	alpha setting for Canny-Deriche filter
     * @param upper	Upper treshold for double tresholding
     * @param lower	lower treshold for double tresholding
     * @return ImageProcessor with 16bit greyscale filtered image
     *
     */
    public static ImageProcessor areaEdge(ImageProcessor ip,
            double radius,
            float alpha,
            float upper,
            float lower) {
        int type = getType(ip);
        if (type == OTHER) {
            IJ.error("area detection", "No action taken. Greyscale or pseudocolored images required!");
            return ip;
        }
        ImageProcessor ipcopy = getDeriche(ip, alpha, radius);
        ipcopy = trin(ipcopy, upper, lower);
        ipcopy = hyst(ipcopy);
        return ipcopy;
    }

    /**
     * Canny-Deriche filtering without previous median filter.
     *
     * @param	ip	ImageProcessor
     * @param	alpha	alpha value for Deriche filter
     * @return	ImageProcessor with 16bit greyscale filtered image.
     *
     */
    public static ImageProcessor getDeriche(ImageProcessor ip, float alpha) {
        return getDeriche(ip, alpha, 0);

    }

    /**
     * Canny-Deriche filtering with previous median filter.
     *
     * @param	ip	ImageProcessor
     * @param	alpha	alpha value for Deriche filter
     * @param radius	radius for the median filter
     * @return	ImageProcessor with 16bit greyscale filtered image.
     *
     */
    public static ImageProcessor getDeriche(ImageProcessor ip, float alpha, double radius) {
        int type = getType(ip);
        if (type == OTHER) {
            IJ.error("Deriche", "No action taken. Greyscale or pseudocolored image required!");
            return ip;
        }
        ArrayList<double[]> arrays = null;
        ImageProcessor ip2 = ip.duplicate();
        RankFilters filter = new RankFilters();
        if (radius > 0) {
            filter.rank(ip2, radius, RankFilters.MEDIAN);
        }

        arrays = dericheCalc(ip2, alpha);
        double[] norm = arrays.get(0);
        double[] angle = arrays.get(1);
        FloatProcessor normfp = new FloatProcessor(ip2.getWidth(), ip2.getHeight(), norm);
        normfp.resetMinAndMax();
        FloatProcessor anglefp = new FloatProcessor(ip2.getWidth(), ip2.getHeight(), angle);
        anglefp.resetMinAndMax();
        ip2 = nonMaximalSuppression(normfp, anglefp);
        return ip2;
    }

    /**
     * Canny-Deriche filtering without previous median filter. Angle output
     *
     * @param	ip	ImageProcessor
     * @param	alpha	alpha value for Deriche filter
     * @return	FloatProcessor with angle image
     *
     */
    public static FloatProcessor getDericheAngle(ImageProcessor ip, float alpha) {
        return getDericheAngle(ip, alpha, 0);
    }

    /**
     * Canny-Deriche filtering with previous median filter.
     *
     * @param	ip	ImageProcessor
     * @param	alpha	alpha value for Deriche filter
     * @param radius	radius for the median filter
     * @return	FloatProcessor with angle image
     *
     */
    public static FloatProcessor getDericheAngle(ImageProcessor ip, float alpha, double radius) {
        int type = getType(ip);
        if (type == OTHER) {
            IJ.error("Deriche", "No action taken. Greyscale or pseudocolored image required");
            return (FloatProcessor) ip;
        }
        ImageProcessor ip2 = ip.duplicate();
        RankFilters filter = new RankFilters();
        if (radius > 0) {
            filter.rank(ip2, radius, RankFilters.MEDIAN);
        }

        double[] angle = dericheCalc(ip2, alpha).get(1);
//        FloatProcessor anglefp = new FloatProcessor(ip2.getWidth(), ip2.getHeight(), angle.getArray()[0]);
//        anglefp.resetMinAndMax();
        return null;//anglefp;
    }

    /**
     * Canny-Deriche filtering without previous median filter. Norm output
     *
     * @param	ip	ImageProcessor
     * @param	alpha	alpha value for Deriche filter
     * @return	FloatProcessor with Norm output
     *
     */
    public static FloatProcessor getDericheNorm(ImageProcessor ip, float alpha) {
        return getDericheNorm(ip, alpha, 0);
    }

    /**
     * Canny-Deriche filtering with previous median filter.
     *
     * @param	ip	ImageProcessor
     * @param	alpha	alpha value for Deriche filter
     * @param radius	radius for the median filter
     * @return	FloatProcessor with Norm output
     *
     */
    public static FloatProcessor getDericheNorm(ImageProcessor ip, float alpha, double radius) {
        int type = getType(ip);
        if (type == OTHER) {
            IJ.error("Deriche", "greyscale or pseudocolored images required");
            return (FloatProcessor) ip;
        }
        ImageProcessor ip2 = ip.duplicate();
        RankFilters filter = new RankFilters();
        if (radius > 0) {
            filter.rank(ip2, radius, RankFilters.MEDIAN);
        }

        double[] norm = dericheCalc(ip2, alpha).get(0);
        FloatProcessor normfp = new FloatProcessor(ip2.getWidth(), ip2.getHeight(), norm);
        normfp.resetMinAndMax();
        return null;//normfp;
    }

    /**
     * Suppression of non local-maxima
     *
     * @param grad the norm gradient image
     * @param ang the angle gradient image
     * @return the image with non local-maxima suppressed
     */
//-------------------------TRESHOLDING METHODS---------------------------------------//	
    /**
     * Double thresholding
     *
     * @param ima original image
     * @param T1 high threshold
     * @param T2 low threshold
     * @return "trinarised" image
     */
    public static ImageProcessor trin(ImageProcessor ima, float T1, float T2) {
        int la = ima.getWidth();
        int ha = ima.getHeight();
        ByteProcessor res = new ByteProcessor(la, ha);
        float pix;

        for (int x = 0; x < la; x++) {
            for (int y = 0; y < ha; y++) {
                pix = ima.getPixelValue(x, y);
                if (pix >= T1) {
                    res.putPixel(x, y, 255);
                } else if (pix >= T2) {
                    res.putPixel(x, y, 128);
                }
            }
        }
        return res;
    }

    /**
     * Hysteresis thresholding
     *
     * @param ima original image
     * @return thresholded image
     */
    public static ImageProcessor hyst(ImageProcessor ima) {
        int la = ima.getWidth();
        int ha = ima.getHeight();
        ImageProcessor res = ima.duplicate();
        boolean change = true;

        // connection
        while (change) {
            change = false;
            for (int x = 1; x < la - 1; x++) {
                for (int y = 1; y < ha - 1; y++) {
                    if (res.getPixelValue(x, y) == 255) {
                        if (res.getPixelValue(x + 1, y) == 128) {
                            change = true;
                            res.putPixelValue(x + 1, y, 255);
                        }
                        if (res.getPixelValue(x - 1, y) == 128) {
                            change = true;
                            res.putPixelValue(x - 1, y, 255);
                        }
                        if (res.getPixelValue(x, y + 1) == 128) {
                            change = true;
                            res.putPixelValue(x, y + 1, 255);
                        }
                        if (res.getPixelValue(x, y - 1) == 128) {
                            change = true;
                            res.putPixelValue(x, y - 1, 255);
                        }
                        if (res.getPixelValue(x + 1, y + 1) == 128) {
                            change = true;
                            res.putPixelValue(x + 1, y + 1, 255);
                        }
                        if (res.getPixelValue(x - 1, y - 1) == 128) {
                            change = true;
                            res.putPixelValue(x - 1, y - 1, 255);
                        }
                        if (res.getPixelValue(x - 1, y + 1) == 128) {
                            change = true;
                            res.putPixelValue(x - 1, y + 1, 255);
                        }
                        if (res.getPixelValue(x + 1, y - 1) == 128) {
                            change = true;
                            res.putPixelValue(x + 1, y - 1, 255);
                        }
                    }
                }
            }
            if (change) {
                for (int x = la - 2; x > 0; x--) {
                    for (int y = ha - 2; y > 0; y--) {
                        if (res.getPixelValue(x, y) == 255) {
                            if (res.getPixelValue(x + 1, y) == 128) {
                                change = true;
                                res.putPixelValue(x + 1, y, 255);
                            }
                            if (res.getPixelValue(x - 1, y) == 128) {
                                change = true;
                                res.putPixelValue(x - 1, y, 255);
                            }
                            if (res.getPixelValue(x, y + 1) == 128) {
                                change = true;
                                res.putPixelValue(x, y + 1, 255);
                            }
                            if (res.getPixelValue(x, y - 1) == 128) {
                                change = true;
                                res.putPixelValue(x, y - 1, 255);
                            }
                            if (res.getPixelValue(x + 1, y + 1) == 128) {
                                change = true;
                                res.putPixelValue(x + 1, y + 1, 255);
                            }
                            if (res.getPixelValue(x - 1, y - 1) == 128) {
                                change = true;
                                res.putPixelValue(x - 1, y - 1, 255);
                            }
                            if (res.getPixelValue(x - 1, y + 1) == 128) {
                                change = true;
                                res.putPixelValue(x - 1, y + 1, 255);
                            }
                            if (res.getPixelValue(x + 1, y - 1) == 128) {
                                change = true;
                                res.putPixelValue(x + 1, y - 1, 255);
                            }
                        }
                    }
                }
            }
        }
        // suppression
        for (int x = 0; x < la; x++) {
            for (int y = 0; y < ha; y++) {
                if (res.getPixelValue(x, y) == 128) {
                    res.putPixelValue(x, y, 0);
                }
            }
        }
        return res;
    }

// ----------------------------PRIVATE METHODS-------------------------------- //
    /**
     * Static method returning an ArrayList of double[] containing the norm
     * array(index 0) and the angle array(index 1).
     *
     */
    private static class OneMatrixThread implements Runnable {

        public OneMatrixThread(int lines, int columns) {
            Matrix a1 = Matrices.random(lines, columns);
        }

        @Override
        public void run() {
        }
    }

    public static ArrayList<double[]> dericheCalc(ImageProcessor ip, float alphaD) {
        double[] norm_deriche = null;
        double[] angle_deriche = null;

        int lines = ip.getHeight();
        int columns = ip.getWidth();
//        float[] nf_grx = null;
        float[] nf_gry = null;
//        int[] a1 = null;
        float[] a2 = new float[lines * columns];
        float[] a3 = new float[lines * columns];
        float[] a4 = new float[lines * columns];

//        int icolonnes = 0;
//        nmem = lines * columns;

        int lineMinus1 = lines - 1;
        int lineMinus2 = lines - 2;
        int lineMinus3 = lines - 3;
        int colMinus1 = columns - 1;
        int colMinus2 = columns - 2;
        int colMinus3 = columns - 3;
//        int icol_1;
//        int icol_2;
        float ad1, ad2, an1, an2, an3, an4, an11;
        /*
         * alloc temporary buffers
         */
        norm_deriche = new double[lines * columns];
        angle_deriche = new double[lines * columns];
        ArrayList<double[]> result = new ArrayList<double[]>();

//        nf_grx = new float[nmem];
        nf_gry = new float[lines * columns];

//        a1 = new int[nmem];

        ad1 = (float) -Math.exp(-alphaD);
        ad2 = 0;
        an1 = 1;
        an2 = 0;
        an3 = (float) Math.exp(-alphaD);
        an4 = 0;
        an11 = 1;

        /*
         * FIRST STEP: Y GRADIENT
         */
        /*
         * x-smoothing
         */
//        for (i = 0; i < lines; i++) {
//            for (j = 0; j < columns; j++) {
//                a1[i * columns + j] = (int) ip.getPixel(j, i);
//
//            }
//        }

        for (int i = 0; i < lines; ++i) {
//            icolonnes = i * columns;
//            icol_1 = i * columns - 1;
//            icol_2 = i * columns - 2;
            a2[i * columns] = an1 * (int) ip.getPixel(0, i);//a1[i*columns+0];
            a2[i * columns + 1] = an1 * (int) ip.getPixel(1, i)//a1[i*columns + 1]
                    + an2 * (int) ip.getPixel(0, i)//a1[icolonnes]
                    - ad1 * a2[i * columns];
            for (int j = 2; j < columns; ++j) {
                a2[i * columns + j] = an1 * (int) ip.getPixel(j, i)//a1[icolonnes + j]
                        + an2 * (int) ip.getPixel(j, i - 1)//a1[icol_1 + j]=a1[i * columns - 1 + j]
                        - ad1 * a2[i * columns - 1 + j] - ad2 * a2[i * columns - 2 + j];
            }
        }

        for (int i = 0; i < lines; ++i) {
//            icolonnes = i * columns;
//            icol_1 = i * columns + 1;
//            icol_2 = i * columns + 2;
            a3[i * columns + colMinus1] = 0;
            a3[i * columns + colMinus2] = an3 * (int) ip.getPixel(colMinus1, i);//a1[icolonnes + col_1]=a1[i * columns + col_1]
            for (int j = colMinus3; j >= 0; --j) {
                a3[i * columns + j] = an3 * (int) ip.getPixel(j, i + 1) //a1[icol_1 + j] 
                        + an4 * (int) ip.getPixel(j, i + 2)//a1[icol_2 + j]
                        - ad1 * a3[i * columns + 1 + j] - ad2 * a3[i * columns + 2 + j];
            }
        }

//        icol_1 = lines * columns;

        for (int i = 0; i < lines * columns; ++i) {
            a2[i] += a3[i];
        }

        /*
         * FIRST STEP Y-GRADIENT : y-derivative
         */
        /*
         * columns top - downn
         */
        for (int j = 0; j < columns; ++j) {
            a3[j] = 0;
            a3[columns + j] = an11 * a2[j] - ad1 * a3[j];
            for (int i = 2; i < lines; ++i) {
                a3[i * columns + j] = an11 * a2[(i - 1) * columns + j]
                        - ad1 * a3[(i - 1) * columns + j] - ad2 * a3[(i - 2) * columns + j];
            }
        }

        /*
         * columns down top
         */
        for (int j = 0; j < columns; ++j) {
            a4[lineMinus1 * columns + j] = 0;
            a4[(lineMinus2 * columns) + j] = -an11 * a2[lineMinus1 * columns + j]
                    - ad1 * a4[lineMinus1 * columns + j];
            for (int i = lineMinus3; i >= 0; --i) {
                a4[i * columns + j] = -an11 * a2[(i + 1) * columns + j]
                        - ad1 * a4[(i + 1) * columns + j] - ad2 * a4[(i + 2) * columns + j];
            }
        }

        for (int i = 0; i < columns * lines; ++i) {
            a3[i] += a4[i];
        }

        for (int i = 0; i < lines; ++i) {
            for (int j = 0; j < columns; ++j) {
                nf_gry[i * columns + j] = a3[i * columns + j];
            }
        }

        /*
         * SECOND STEP X-GRADIENT
         */
//        for (int i = 0; i < lines; ++i) {
//            for (int j = 0; j < columns; ++j) {
//                a1[i * columns + j] = (int) (ip.getPixel(j, i));
//            }
//        }

        for (int i = 0; i < lines; ++i) {
//            icolonnes = i * columns;
//            icol_1 = i * columns - 1;
//            icol_2 = i * columns - 2;
            a2[i * columns] = 0;
            a2[i * columns + 1] = an11 * (int) ip.getPixel(0, i);//a1[icolonnes];
            for (int j = 2; j < columns; ++j) {
                a2[i * columns + j] = an11 * (int) ip.getPixel(j, i - 1)//a1[icol_1 + j]
                        - ad1 * a2[i * columns - 1 + j] - ad2 * a2[i * columns - 2 + j];
            }
        }

        for (int i = 0; i < lines; ++i) {
//            icolonnes = i * columns;
//            icol_1 = i * columns + 1;
//            icol_2 = i * columns + 2;
            a3[i * columns + colMinus1] = 0;
            a3[i * columns + colMinus2] = -an11 * (int) ip.getPixel(colMinus1, i);//a1[icolonnes + col_1];
            for (int j = colMinus3; j >= 0; --j) {
                a3[i * columns + j] = -an11 * (int) ip.getPixel(j, i + 1)//a1[icol_1 + j]
                        - ad1 * a3[i * columns - 1 + j] - ad2 * a3[i * columns - 2 + j];
            }
        }
//        icol_1 = lines * columns;
        for (int i = 0; i < lines * columns; ++i) {
            a2[i] += a3[i];
        }

        /*
         * on the columns
         */
        /*
         * columns top down
         */
        for (int j = 0; j < columns; ++j) {
            a3[j] = an1 * a2[j];
            a3[columns + j] = an1 * a2[columns + j] + an2 * a2[j]
                    - ad1 * a3[j];
            for (int i = 2; i < lines; ++i) {
                a3[i * columns + j] = an1 * a2[i * columns + j] + an2 * a2[(i - 1) * columns + j]
                        - ad1 * a3[(i - 1) * columns + j] - ad2 * a3[(i - 2) * columns + j];
            }
        }

        /*
         * columns down top
         */
        for (int j = 0; j < columns; ++j) {
            a4[lineMinus1 * columns + j] = 0;
            a4[lineMinus2 * columns + j] = an3 * a2[lineMinus1 * columns + j] - ad1 * a4[lineMinus1 * columns + j];
            for (int i = lineMinus3; i >= 0; --i) {
                a4[i * columns + j] = an3 * a2[(i + 1) * columns + j] + an4 * a2[(i + 2) * columns + j]
                        - ad1 * a4[(i + 1) * columns + j] - ad2 * a4[(i + 2) * columns + j];
            }
        }

        for (int i = 0; i < columns * lines; ++i) {
            a3[i] += a4[i];
        }

//        for (int i = 0; i < lines; i++) {
//            for (int j = 0; j < columns; j++) {
//                nf_grx[i * columns + j] = a3[i * columns + j];
//            }
//        }
        for (int i = 0; i < columns * lines; ++i) {
            norm_deriche[i] = modul(a3[i], nf_gry[i]);
            angle_deriche[i] = angle(a3[i], nf_gry[i]);
        }
        /*
         * SECOND STEP X-GRADIENT : the x-gradient is done
         */
        /*
         * THIRD STEP : NORM
         */
        /*
         * computation of the magnitude and angle
         */
//        for (int i = 0; i < lines; i++) {
//            for (int j = 0; j < columns; j++) {
//                a2[i * columns + j] = nf_gry[i * columns + j];
//            }
//        }
//        icol_1 = columns * lines;


        result.add(norm_deriche);
        result.add(angle_deriche);
        return result;
    }

    /**
     * Suppression of non local-maxima
     *
     * @param grad the norm gradient image
     * @param ang the angle gradient image
     * @return the image with non local-maxima suppressed
     */
    public static ImageProcessor nonMaximalSuppression(ImageProcessor grad, ImageProcessor ang) {
        FloatProcessor res = new FloatProcessor(grad.getWidth(), grad.getHeight());

        int la = grad.getWidth();
        int ha = grad.getHeight();

        float ag;
        float pix1 = 0;
        float pix2 = 0;
        float pix;

        for (int x = 1; x < la - 1; x++) {
            for (int y = 1; y < ha - 1; y++) {
                ag = ang.getPixelValue(x, y);
                if ((ag > -22.5) && (ag <= 22.5)) {
                    pix1 = grad.getPixelValue(x + 1, y);
                    pix2 = grad.getPixelValue(x - 1, y);
                } else if ((ag > 22.5) && (ag <= 67.5)) {
                    pix1 = grad.getPixelValue(x + 1, y - 1);
                    pix2 = grad.getPixelValue(x - 1, y + 1);
                } else if (((ag > 67.5) && (ag <= 90)) || ((ag < -67.5) && (ag >= -90))) {
                    pix1 = grad.getPixelValue(x, y - 1);
                    pix2 = grad.getPixelValue(x, y + 1);
                } else if ((ag < -22.5) && (ag >= -67.5)) {
                    pix1 = grad.getPixelValue(x + 1, y + 1);
                    pix2 = grad.getPixelValue(x - 1, y - 1);
                }
                pix = grad.getPixelValue(x, y);
                if ((pix >= pix1) && (pix >= pix2)) {
                    res.putPixelValue(x, y, pix);
                }
            }
        }
        return res;
    }

    /**
     * modul
     *
     * @param dx derivative in x
     * @param dy derivative in y
     * @return norm of gradient
     */
    public static double modul(float dx, float dy) {
        return (Math.sqrt(dx * dx + dy * dy));
    }

    /**
     * angle
     *
     * @param dx derivative in x
     * @param dy derivative in y
     * @return angle of gradient
     */
    public static double angle(float dx, float dy) {
        return (-Math.toDegrees(Math.atan(dy / dx)));
    }

    /**
     * Checks the type of the image
     *
     */
    static int getType(ImageProcessor ip) {
        if (!typed) {
            typed = true;
            if (ip instanceof ByteProcessor & !ip.isColorLut()) {
                type = BYTE;
            } else if (ip instanceof ShortProcessor) {
                type = SHORT;
            } else if (ip instanceof FloatProcessor) {
                type = FLOAT;
            } else {
                type = OTHER;
            }
        }
        return type;
    }
}
