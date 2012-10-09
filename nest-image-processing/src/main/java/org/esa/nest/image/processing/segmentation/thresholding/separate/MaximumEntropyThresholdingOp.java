/*
 * Copyright (C) 2012 by Array Systems Computing Inc. http://www.array.ca
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
package org.esa.nest.image.processing.segmentation.thresholding.separate;

import com.bc.ceres.core.ProgressMonitor;
import ij.ImagePlus;
import ij.process.AutoThresholder;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;
import java.awt.Rectangle;
import java.awt.image.BufferedImage;
import java.awt.image.RenderedImage;
import java.util.HashMap;
import java.util.Map;
import javax.swing.JOptionPane;
import org.esa.beam.framework.datamodel.Band;
import org.esa.beam.framework.datamodel.Product;
import org.esa.beam.framework.datamodel.ProductData;
import org.esa.beam.framework.gpf.Operator;
import org.esa.beam.framework.gpf.OperatorException;
import org.esa.beam.framework.gpf.OperatorSpi;
import org.esa.beam.framework.gpf.Tile;
import org.esa.beam.framework.gpf.annotations.OperatorMetadata;
import org.esa.beam.framework.gpf.annotations.Parameter;
import org.esa.beam.framework.gpf.annotations.SourceProduct;
import org.esa.beam.framework.gpf.annotations.TargetProduct;
import org.esa.nest.gpf.OperatorUtils;

@OperatorMetadata(alias = "MaximumEntropyThresholding",
category = "SAR Tools\\Image Processing",
description = "MaximumEntropyThresholding")
public class MaximumEntropyThresholdingOp extends Operator {

    public static float[] probabilityHistogram;
    final static int MAX_VALUE = 256;
    final static int MIN_VALUE = 0;
    public static int N;
    @SourceProduct(alias = "source")
    private Product sourceProduct = null;
    @TargetProduct
    private Product targetProduct;
    @Parameter(description = "The list of source bands.", alias = "sourceBands", itemAlias = "band",
    rasterDataNodeType = Band.class, label = "Source Bands")
    private String[] sourceBandNames;
    @Parameter(description = "HighThreshold", defaultValue = "100", label = "HighThreshold")
    private float highThreshold = 100f;
    @Parameter(description = "LowThreshold", defaultValue = "10", label = "LowThreshold")
    private float lowThreshold = 10f;
    private final Map<String, String[]> targetBandNameToSourceBandName = new HashMap<String, String[]>();
    private int sourceImageWidth;
    private int sourceImageHeight;
    private boolean processed = false;
    private int halfSizeX;
    private int halfSizeY;
    private int filterSizeX = 3;
    private int filterSizeY = 3;
    private static ImagePlus fullImagePlus;
    private static ByteProcessor fullByteProcessor;

    /**
     * Initializes this operator and sets the one and only target product.
     * <p>The target product can be either defined by a field of type {@link org.esa.beam.framework.datamodel.Product}
     * annotated with the
     * {@link org.esa.beam.framework.gpf.annotations.TargetProduct TargetProduct}
     * annotation or by calling {@link #setTargetProduct} method.</p> <p>The
     * framework calls this method after it has created this operator. Any
     * client code that must be performed before computation of tile data should
     * be placed here.</p>
     *
     * @throws org.esa.beam.framework.gpf.OperatorException If an error occurs
     * during operator initialization.
     * @see #getTargetProduct()
     */
    @Override
    public void initialize() throws OperatorException {

        try {
            sourceImageWidth = sourceProduct.getSceneRasterWidth();
            sourceImageHeight = sourceProduct.getSceneRasterHeight();

            halfSizeX = filterSizeX / 2;
            halfSizeY = filterSizeY / 2;

            createTargetProduct();

        } catch (Throwable e) {
            OperatorUtils.catchOperatorException(getId(), e);
        }
    }

    /**
     * Create target product.
     *
     * @throws Exception The exception.
     */
    private void createTargetProduct() throws Exception {

        targetProduct = new Product(sourceProduct.getName(),
                sourceProduct.getProductType(),
                sourceImageWidth,
                sourceImageHeight);

        OperatorUtils.copyProductNodes(sourceProduct, targetProduct);

        OperatorUtils.addSelectedBands(
                sourceProduct, sourceBandNames, targetProduct, targetBandNameToSourceBandName, true);
    }

    /**
     * Called by the framework in order to compute a tile for the given target
     * band. <p>The default implementation throws a runtime exception with the
     * message "not implemented".</p>
     *
     * @param targetBand The target band.
     * @param targetTile The current tile associated with the target band to be
     * computed.
     * @param pm A progress monitor which should be used to determine
     * computation cancelation requests.
     * @throws org.esa.beam.framework.gpf.OperatorException If an error occurs
     * during computation of the target raster.
     */
    @Override
    public void computeTile(Band targetBand, Tile targetTile, ProgressMonitor pm) throws OperatorException {

        try {
            final Rectangle targetTileRectangle = targetTile.getRectangle();
            final int x0 = targetTileRectangle.x;
            final int y0 = targetTileRectangle.y;
            final int w = targetTileRectangle.width;
            final int h = targetTileRectangle.height;

            final Rectangle sourceTileRectangle = getSourceTileRectangle(x0, y0, w, h);
            Tile sourceRaster;
            final String[] srcBandNames = targetBandNameToSourceBandName.get(targetBand.getName());
            Band sourceBand = sourceProduct.getBand(srcBandNames[0]);
            sourceRaster = getSourceTile(sourceBand, sourceTileRectangle);
            if (sourceRaster == null) {
                throw new OperatorException("Cannot get source tile");
            }

            computeMaximumEntropyThresholding(sourceBand, sourceRaster, targetTile, x0, y0, w, h, pm);

        } catch (Throwable e) {
            OperatorUtils.catchOperatorException(getId(), e);
        } finally {
            pm.done();
        }
    }

    /**
     * Apply HysteresisThesholding
     *
     * @param sourceRaster The source tile for the band.
     * @param targetTile The current tile associated with the target band to be
     * computed.
     * @param x0 X coordinate for the upper-left point of the
     * target_Tile_Rectangle.
     * @param y0 Y coordinate for the upper-left point of the
     * target_Tile_Rectangle.
     * @param w Width for the target_Tile_Rectangle.
     * @param h Hight for the target_Tile_Rectangle.
     * @param pm A progress monitor which should be used to determine
     * computation cancellation requests.
     * @throws org.esa.beam.framework.gpf.OperatorException If an error occurs
     * during computation of the filtered value.
     */
    private synchronized void computeMaximumEntropyThresholding(final Band sourceBand, final Tile sourceRaster,
            final Tile targetTile, final int x0, final int y0, final int w, final int h,
            final ProgressMonitor pm) {

        int threshold = 0;
        if (!processed) {
 

            processed = true;
        }

        final Rectangle srcTileRectangle = sourceRaster.getRectangle();

        ImageProcessor aPartProcessor = fullByteProcessor.duplicate();

        aPartProcessor.setRoi(srcTileRectangle);

        ImageProcessor roiImageProcessor = aPartProcessor.crop();

        final ProductData trgData = targetTile.getDataBuffer();
        final ProductData sourceData = ProductData.createInstance((byte[]) roiImageProcessor.getPixels());

        final int maxY = y0 + h;
        final int maxX = x0 + w;
        for (int y = y0; y < maxY; ++y) {
            for (int x = x0; x < maxX; ++x) {

                trgData.setElemFloatAt(targetTile.getDataBufferIndex(x, y),
                        sourceData.getElemFloatAt(sourceRaster.getDataBufferIndex(x, y)));
            }
        }
    }

   

    /**
     * Hysteresis thresholding
     *
     * @param imageProcessor original image
     * @return thresholded image
     */
    ImageProcessor maximumEntropyThresholding(ByteProcessor imageProcessor) {

        ImageProcessor returnedProcessor = imageProcessor.duplicate();

        int[] hist = imageProcessor.getHistogram();
        int threshold = entropySplit(hist);
        returnedProcessor.threshold(threshold);
        return returnedProcessor;
    }

    /**
     * Calculate maximum entropy split of a histogram.
     *
     * @param hist histogram to be thresholded.
     *
     * @return index of the maximum entropy split.`
     */
    private int entropySplit(int[] data) {

        // Implements Kapur-Sahoo-Wong (Maximum Entropy) thresholding method
        // Kapur J.N., Sahoo P.K., and Wong A.K.C. (1985) "A New Method for
        // Gray-Level Picture Thresholding Using the Entropy of the Histogram"
        // Graphical Models and Image Processing, 29(3): 273-285
        // M. Emre Celebi
        // 06.15.2007
        // Ported to ImageJ plugin by G.Landini from E Celebi's fourier_0.8 routines
        int threshold = -1;
        int ih, it;
        int first_bin;
        int last_bin;
        double tot_ent;  /*
         * total entropy
         */
        double max_ent;  /*
         * max entropy
         */
        double ent_back; /*
         * entropy of the background pixels at a given threshold
         */
        double ent_obj;  /*
         * entropy of the object pixels at a given threshold
         */
        double[] norm_histo = new double[data.length]; /*
         * normalized histogram
         */
        double[] P1 = new double[data.length]; /*
         * cumulative normalized histogram
         */
        double[] P2 = new double[data.length];

        int total = 0;
        for (ih = 0; ih < data.length; ih++) {
            total += data[ih];
        }

        for (ih = 0; ih < data.length; ih++) {
            norm_histo[ih] = (double) data[ih] / total;
        }

        P1[0] = norm_histo[0];
        P2[0] = 1.0 - P1[0];
        for (ih = 1; ih < data.length; ih++) {
            P1[ih] = P1[ih - 1] + norm_histo[ih];
            P2[ih] = 1.0 - P1[ih];
        }

        /*
         * Determine the first non-zero bin
         */
        first_bin = 0;
        for (ih = 0; ih < data.length; ih++) {
            if (!(Math.abs(P1[ih]) < 2.220446049250313E-16)) {
                first_bin = ih;
                break;
            }
        }

        /*
         * Determine the last non-zero bin
         */
        last_bin = data.length - 1;
        for (ih = data.length - 1; ih >= first_bin; ih--) {
            if (!(Math.abs(P2[ih]) < 2.220446049250313E-16)) {
                last_bin = ih;
                break;
            }
        }

        // Calculate the total entropy each gray-level
        // and find the threshold that maximizes it 
        max_ent = Double.MIN_VALUE;

        for (it = first_bin; it <= last_bin; it++) {
            /*
             * Entropy of the background pixels
             */
            ent_back = 0.0;
            for (ih = 0; ih <= it; ih++) {
                if (data[ih] != 0) {
                    ent_back -= (norm_histo[ih] / P1[it]) * Math.log(norm_histo[ih] / P1[it]);
                }
            }

            /*
             * Entropy of the object pixels
             */
            ent_obj = 0.0;
            for (ih = it + 1; ih < data.length; ih++) {
                if (data[ih] != 0) {
                    ent_obj -= (norm_histo[ih] / P2[it]) * Math.log(norm_histo[ih] / P2[it]);
                }
            }

            /*
             * Total entropy
             */
            tot_ent = ent_back + ent_obj;

            // IJ.log(""+max_ent+"  "+tot_ent);
            if (max_ent < tot_ent) {
                max_ent = tot_ent;
                threshold = it;
            }
        }
        return threshold;
    }

    /**
     * Get source tile rectangle.
     *
     * @param x0 X coordinate of the upper left corner point of the target tile
     * rectangle.
     * @param y0 Y coordinate of the upper left corner point of the target tile
     * rectangle.
     * @param w The width of the target tile rectangle.
     * @param h The height of the target tile rectangle.
     * @return The source tile rectangle.
     */
    private Rectangle getSourceTileRectangle(int x0, int y0, int w, int h) {

        int sx0 = x0;
        int sy0 = y0;
        int sw = w;
        int sh = h;

        if (x0 >= halfSizeX) {
            sx0 -= halfSizeX;
            sw += halfSizeX;
        }

        if (y0 >= halfSizeY) {
            sy0 -= halfSizeY;
            sh += halfSizeY;
        }

        if (x0 + w + halfSizeX <= sourceImageWidth) {
            sw += halfSizeX;
        }

        if (y0 + h + halfSizeY <= sourceImageHeight) {
            sh += halfSizeY;
        }

        return new Rectangle(sx0, sy0, sw, sh);
    }

    @Override
    public void dispose() {
        fullByteProcessor = null;
        fullImagePlus = null;
        processed = false;
    }

    /**
     * The SPI is used to register this operator in the graph processing
     * framework via the SPI configuration file
     * {@code META-INF/services/org.esa.beam.framework.gpf.OperatorSpi}. This
     * class may also serve as a factory for new operator instances.
     *
     * @see OperatorSpi#createOperator()
     * @see OperatorSpi#createOperator(java.util.Map, java.util.Map)
     */
    public static class Spi extends OperatorSpi {

        public Spi() {
            super(MaximumEntropyThresholdingOp.class);
            setOperatorUI(MaximumEntropyThresholdingOpUI.class);
        }
    }
}
