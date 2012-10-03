package org.esa.nest.image.processing.utils.image;

import java.util.ArrayList;

import ij.*;
import ij.process.*;
import ij.plugin.filter.RankFilters;

/**
 * Class containing some static methods for Edge Detection, using median filter,
 * Canny-Deriche filtering, double thresholding, hysteresis thresholding. These
 * are meant to be implemented in another plug-in, hence the basic
 * functionality. Methods based on the plug-ins deriche_ and hysteresis_
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
     * @param imageProcessor	ImageProcessor
     * @param radius	radius for the median filter
     * @param alpha	alpha setting for Canny-Deriche filter
     * @param upper	Upper threshold for double thresholding
     * @param lower	lower threshold for double thresholding
     * @return ImageProcessor with 16bit greyscale filtered image
     *
     */
    public static ImageProcessor areaEdge(ImageProcessor imageProcessor,
            double radius,
            float alpha,
            float upper,
            float lower) {
        int ipType = getType(imageProcessor);
        if (ipType == OTHER) {
            IJ.error("area detection", "No action taken. Greyscale or pseudocolored images required!");
            return imageProcessor;
        }
        ImageProcessor ipcopy = getDeriche(imageProcessor, alpha, radius);
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
    public static ImageProcessor getDeriche(ImageProcessor imageProcessor, float alpha) {
        return getDeriche(imageProcessor, alpha, 0);

    }

    /**
     * Canny-Deriche filtering with previous median filter.
     *
     * @param	imageProcessor	ImageProcessor
     * @param	alpha	alpha value for Deriche filter
     * @param radius	radius for the median filter
     * @return	ImageProcessor with 16bit greyscale filtered image.
     *
     */
    public static ImageProcessor getDeriche(ImageProcessor imageProcessor, float alpha, double radius) {
        int ipType = getType(imageProcessor);
        if (ipType == OTHER) {
            IJ.error("Deriche", "No action taken. Greyscale or pseudocolored image required!");
            return imageProcessor;
        }
        ArrayList<double[]> arrays = null;
        ImageProcessor ip2 = imageProcessor.duplicate();
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
    public static FloatProcessor getDericheAngle(ImageProcessor imageProcessor, float alpha) {
        return getDericheAngle(imageProcessor, alpha, 0);
    }

    /**
     * Canny-Deriche filtering with previous median filter.
     *
     * @param	imageProcessor	ImageProcessor
     * @param	alpha	alpha value for Deriche filter
     * @param radius	radius for the median filter
     * @return	FloatProcessor with angle image
     *
     */
    public static FloatProcessor getDericheAngle(ImageProcessor imageProcessor, float alpha, double radius) {
        int ipType = getType(imageProcessor);
        if (ipType == OTHER) {
            IJ.error("Deriche", "No action taken. Greyscale or pseudocolored image required");
            return (FloatProcessor) imageProcessor;
        }
        ImageProcessor ip2 = imageProcessor.duplicate();
        RankFilters filter = new RankFilters();
        if (radius > 0) {
            filter.rank(ip2, radius, RankFilters.MEDIAN);
        }

        double[] angle = dericheCalc(ip2, alpha).get(1);
        FloatProcessor anglefp = new FloatProcessor(ip2.getWidth(), ip2.getHeight(), angle);
        anglefp.resetMinAndMax();
        return anglefp;
    }

    /**
     * Canny-Deriche filtering without previous median filter. Norm output
     *
     * @param	imageProcessor	ImageProcessor
     * @param	alpha	alpha value for Deriche filter
     * @return	FloatProcessor with Norm output
     *
     */
    public static FloatProcessor getDericheNorm(ImageProcessor imageProcessor, float alpha) {
        return getDericheNorm(imageProcessor, alpha, 0);
    }

    /**
     * Canny-Deriche filtering with previous median filter.
     *
     * @param	imageProcessor	ImageProcessor
     * @param	alpha	alpha value for Deriche filter
     * @param radius	radius for the median filter
     * @return	FloatProcessor with Norm output
     *
     */
    public static FloatProcessor getDericheNorm(ImageProcessor imageProcessor, float alpha, double radius) {
        int ipType = getType(imageProcessor);
        if (ipType == OTHER) {
            IJ.error("Deriche", "greyscale or pseudocolored images required");
            return (FloatProcessor) imageProcessor;
        }
        ImageProcessor ip2 = imageProcessor.duplicate();
        RankFilters filter = new RankFilters();
        if (radius > 0) {
            filter.rank(ip2, radius, RankFilters.MEDIAN);
        }

        double[] norm = dericheCalc(ip2, alpha).get(0);
        FloatProcessor normfp = new FloatProcessor(ip2.getWidth(), ip2.getHeight(), norm);
        normfp.resetMinAndMax();
        return normfp;
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
    public static ImageProcessor trin(ImageProcessor imageProcessor, float T1, float T2) {
        int width = imageProcessor.getWidth();
        int height = imageProcessor.getHeight();
        ByteProcessor resultProcessor = new ByteProcessor(width, height);
        float pix;

        for (int x = 0; x < width; x++) {
            for (int y = 0; y < height; y++) {
                pix = imageProcessor.getPixelValue(x, y);
                if (pix >= T1) {
                    resultProcessor.putPixel(x, y, 255);
                } else if (pix >= T2) {
                    resultProcessor.putPixel(x, y, 128);
                }
            }
        }
        return resultProcessor;
    }

    /**
     * Hysteresis thresholding
     *
     * @param ima original image
     * @return thresholded image
     */
    public static ImageProcessor hyst(ImageProcessor imageProcessor) {
        int la = imageProcessor.getWidth();
        int ha = imageProcessor.getHeight();
        ImageProcessor res = imageProcessor.duplicate();
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
    public static ArrayList<double[]> dericheCalc(ImageProcessor imageProcessor, float alphaD) {
        double[] norm_deriche = null;
        double[] angle_deriche = null;

        int nmem;
        float[] nf_grx = null;
        float[] nf_gry = null;
        int[] a1 = null;
        float[] a2 = null;
        float[] a3 = null;
        float[] a4 = null;

        int iHeight = 0;
        int height;
        int width;
        int lig_1;
        int lig_2;
        int lig_3;
        int col_1;
        int col_2;
        int col_3;
        int icol_1;
        int icol_2;
        int i;
        int j;
        float ad1;
        float ad2;
        float an1;
        float an2;
        float an3;
        float an4;
        float an11;

        height = imageProcessor.getHeight();
        width = imageProcessor.getWidth();
        nmem = height * width;

        lig_1 = height - 1;
        lig_2 = height - 2;
        lig_3 = height - 3;
        col_1 = width - 1;
        col_2 = width - 2;
        col_3 = width - 3;

        /*
         * alloc temporary buffers
         */
        norm_deriche = new double[nmem];
        angle_deriche = new double[nmem];
        ArrayList<double[]> result = new ArrayList<double[]>();

        nf_grx = new float[nmem];
        nf_gry = new float[nmem];

        a1 = new int[nmem];
        a2 = new float[nmem];
        a3 = new float[nmem];
        a4 = new float[nmem];

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
        for (i = 0; i < height; i++) {
            for (j = 0; j < width; j++) {
                a1[i * width + j] = (int) imageProcessor.getPixel(j, i);

            }
        }

        for (i = 0; i < height; ++i) {
            iHeight = i * width;
            icol_1 = iHeight - 1;
            icol_2 = iHeight - 2;
            a2[iHeight] = an1 * a1[iHeight];
            a2[iHeight + 1] = an1 * a1[iHeight + 1]
                    + an2 * a1[iHeight] - ad1 * a2[iHeight];
            for (j = 2; j < width; ++j) {
                a2[iHeight + j] = an1 * a1[iHeight + j] + an2 * a1[icol_1 + j]
                        - ad1 * a2[icol_1 + j] - ad2 * a2[icol_2 + j];
            }
        }

        for (i = 0; i < height; ++i) {
            iHeight = i * width;
            icol_1 = iHeight + 1;
            icol_2 = iHeight + 2;
            a3[iHeight + col_1] = 0;
            a3[iHeight + col_2] = an3 * a1[iHeight + col_1];
            for (j = col_3; j >= 0; --j) {
                a3[iHeight + j] = an3 * a1[icol_1 + j] + an4 * a1[icol_2 + j]
                        - ad1 * a3[icol_1 + j] - ad2 * a3[icol_2 + j];
            }
        }

        icol_1 = height * width;

        for (i = 0; i < icol_1; ++i) {
            a2[i] += a3[i];
        }

        /*
         * FIRST STEP Y-GRADIENT : y-derivative
         */
        /*
         * columns top - downn
         */
        for (j = 0; j < width; ++j) {
            a3[j] = 0;
            a3[width + j] = an11 * a2[j] - ad1 * a3[j];
            for (i = 2; i < height; ++i) {
                a3[i * width + j] = an11 * a2[(i - 1) * width + j]
                        - ad1 * a3[(i - 1) * width + j] - ad2 * a3[(i - 2) * width + j];
            }
        }

        /*
         * columns down top
         */
        for (j = 0; j < width; ++j) {
            a4[lig_1 * width + j] = 0;
            a4[(lig_2 * width) + j] = -an11 * a2[lig_1 * width + j]
                    - ad1 * a4[lig_1 * width + j];
            for (i = lig_3; i >= 0; --i) {
                a4[i * width + j] = -an11 * a2[(i + 1) * width + j]
                        - ad1 * a4[(i + 1) * width + j] - ad2 * a4[(i + 2) * width + j];
            }
        }

        icol_1 = width * height;
        for (i = 0; i < icol_1; ++i) {
            a3[i] += a4[i];
        }

        for (i = 0; i < height; ++i) {
            for (j = 0; j < width; ++j) {
                nf_gry[i * width + j] = a3[i * width + j];
            }
        }

        /*
         * SECOND STEP X-GRADIENT
         */
        for (i = 0; i < height; ++i) {
            for (j = 0; j < width; ++j) {
                a1[i * width + j] = (int) (imageProcessor.getPixel(j, i));
            }
        }

        for (i = 0; i < height; ++i) {
            iHeight = i * width;
            icol_1 = iHeight - 1;
            icol_2 = iHeight - 2;
            a2[iHeight] = 0;
            a2[iHeight + 1] = an11 * a1[iHeight];
            for (j = 2; j < width; ++j) {
                a2[iHeight + j] = an11 * a1[icol_1 + j]
                        - ad1 * a2[icol_1 + j] - ad2 * a2[icol_2 + j];
            }
        }

        for (i = 0; i < height; ++i) {
            iHeight = i * width;
            icol_1 = iHeight + 1;
            icol_2 = iHeight + 2;
            a3[iHeight + col_1] = 0;
            a3[iHeight + col_2] = -an11 * a1[iHeight + col_1];
            for (j = col_3; j >= 0; --j) {
                a3[iHeight + j] = -an11 * a1[icol_1 + j]
                        - ad1 * a3[icol_1 + j] - ad2 * a3[icol_2 + j];
            }
        }
        icol_1 = height * width;
        for (i = 0; i < icol_1; ++i) {
            a2[i] += a3[i];
        }

        /*
         * on the columns
         */
        /*
         * columns top down
         */
        for (j = 0; j < width; ++j) {
            a3[j] = an1 * a2[j];
            a3[width + j] = an1 * a2[width + j] + an2 * a2[j]
                    - ad1 * a3[j];
            for (i = 2; i < height; ++i) {
                a3[i * width + j] = an1 * a2[i * width + j] + an2 * a2[(i - 1) * width + j]
                        - ad1 * a3[(i - 1) * width + j] - ad2 * a3[(i - 2) * width + j];
            }
        }

        /*
         * columns down top
         */
        for (j = 0; j < width; ++j) {
            a4[lig_1 * width + j] = 0;
            a4[lig_2 * width + j] = an3 * a2[lig_1 * width + j] - ad1 * a4[lig_1 * width + j];
            for (i = lig_3; i >= 0; --i) {
                a4[i * width + j] = an3 * a2[(i + 1) * width + j] + an4 * a2[(i + 2) * width + j]
                        - ad1 * a4[(i + 1) * width + j] - ad2 * a4[(i + 2) * width + j];
            }
        }

        icol_1 = width * height;
        for (i = 0; i < icol_1; ++i) {
            a3[i] += a4[i];
        }

        for (i = 0; i < height; i++) {
            for (j = 0; j < width; j++) {
                nf_grx[i * width + j] = a3[i * width + j];
            }
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
        for (i = 0; i < height; i++) {
            for (j = 0; j < width; j++) {
                a2[i * width + j] = nf_gry[i * width + j];
            }
        }
        icol_1 = width * height;
        for (i = 0; i < icol_1; ++i) {
            norm_deriche[i] = module(nf_grx[i], nf_gry[i]);
            angle_deriche[i] = angle(nf_grx[i], nf_gry[i]);
        }
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
    public static ImageProcessor nonMaximalSuppression(ImageProcessor imageProcessor, ImageProcessor ang) {
        FloatProcessor res = new FloatProcessor(imageProcessor.getWidth(), imageProcessor.getHeight());

        int width = imageProcessor.getWidth();
        int height = imageProcessor.getHeight();

        float ag;
        float pix1 = 0;
        float pix2 = 0;
        float pix;

        for (int x = 1; x < width - 1; x++) {
            for (int y = 1; y < height - 1; y++) {
                ag = ang.getPixelValue(x, y);
                if ((ag > -22.5) && (ag <= 22.5)) {
                    pix1 = imageProcessor.getPixelValue(x + 1, y);
                    pix2 = imageProcessor.getPixelValue(x - 1, y);
                } else if ((ag > 22.5) && (ag <= 67.5)) {
                    pix1 = imageProcessor.getPixelValue(x + 1, y - 1);
                    pix2 = imageProcessor.getPixelValue(x - 1, y + 1);
                } else if (((ag > 67.5) && (ag <= 90)) || ((ag < -67.5) && (ag >= -90))) {
                    pix1 = imageProcessor.getPixelValue(x, y - 1);
                    pix2 = imageProcessor.getPixelValue(x, y + 1);
                } else if ((ag < -22.5) && (ag >= -67.5)) {
                    pix1 = imageProcessor.getPixelValue(x + 1, y + 1);
                    pix2 = imageProcessor.getPixelValue(x - 1, y - 1);
                }
                pix = imageProcessor.getPixelValue(x, y);
                if ((pix >= pix1) && (pix >= pix2)) {
                    res.putPixelValue(x, y, pix);
                }
            }
        }
        return res;
    }

    /**
     * module
     *
     * @param dx derivative in x
     * @param dy derivative in y
     * @return norm of gradient
     */
    public static double module(float dx, float dy) {
        return (Math.sqrt(dx * dx + dy * dy));
    }

    /**
     * angle
     *
     * @param derivativeX derivative in x
     * @param derivativeY derivative in y
     * @return angle of gradient
     */
    public static double angle(float derivativeX, float derivativeY) {
        return (-Math.toDegrees(Math.atan(derivativeY / derivativeX)));
    }

    /**
     * Checks the type of the image
     *
     */
    static int getType(ImageProcessor imageProcessor) {
        if (!typed) {
            typed = true;
            if (imageProcessor instanceof ByteProcessor & !imageProcessor.isColorLut()) {
                type = BYTE;
            } else if (imageProcessor instanceof ShortProcessor) {
                type = SHORT;
            } else if (imageProcessor instanceof FloatProcessor) {
                type = FLOAT;
            } else {
                type = OTHER;
            }
        }
        return type;
    }
}
