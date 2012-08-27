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
package org.esa.nest.image.processing.segmentation;

import com.bc.ceres.core.ProgressMonitor;
import ij.ImagePlus;
import ij.process.ImageProcessor;
import ij.process.ShortProcessor;
import java.awt.Graphics2D;
import java.awt.Image;
import java.awt.Rectangle;
import java.awt.image.BufferedImage;
import java.awt.image.RenderedImage;
import java.awt.image.VolatileImage;
import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.imageio.ImageIO;
import javax.media.jai.PlanarImage;
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

@OperatorMetadata(alias = "OtsuThresholding",
category = "SAR Tools\\Image Processing",
description = "OtsuThresholdingOp")
public class OtsuThresholdingOp extends Operator {

    private static float[] probabilityHistogram;
    public static boolean probabilityHistogramDone;
    public static int N;
    @SourceProduct(alias = "source")
    private Product sourceProduct = null;
    @TargetProduct
    private Product targetProduct;
    @Parameter(description = "The list of source bands.", alias = "sourceBands", itemAlias = "band",
    rasterDataNodeType = Band.class, label = "Source Bands")
    private String[] sourceBandNames;
    private final Map<String, String[]> targetBandNameToSourceBandName = new HashMap<String, String[]>();
    private int sourceImageWidth;
    private int sourceImageHeight;
    private boolean processed = false;

    @Override
    public void initialize() throws OperatorException {
        try {
            sourceImageWidth = sourceProduct.getSceneRasterWidth();
            sourceImageHeight = sourceProduct.getSceneRasterHeight();

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
                sourceProduct, sourceBandNames, targetProduct,
                targetBandNameToSourceBandName, true);
    }

    @Override
    public void computeTile(Band targetBand, Tile targetTile, ProgressMonitor progressMonitor) throws OperatorException {

        try {
            final Rectangle targetTileRectangle = targetTile.getRectangle();
            final int x0 = targetTileRectangle.x;
            final int y0 = targetTileRectangle.y;
            final int w = targetTileRectangle.width;
            final int h = targetTileRectangle.height;

            Tile sourceRaster;
            final String[] srcBandNames = targetBandNameToSourceBandName.get(targetBand.getName());
            Band sourceBand = sourceProduct.getBand(srcBandNames[0]);
            sourceRaster = getSourceTile(sourceBand, targetTileRectangle);

            if (sourceRaster == null) {
                throw new OperatorException("Cannot get source tile");
            }
            if (!processed) {
                compute(sourceBand, progressMonitor);
            }
            final ProductData masterData = sourceRaster.getDataBuffer();
            final ProductData targetData = targetTile.getDataBuffer();
            for (int y = targetTileRectangle.y; y < targetTileRectangle.y + targetTileRectangle.height; y++) {
                for (int x = targetTileRectangle.x; x < targetTileRectangle.x + targetTileRectangle.width; x++) {
                    final int index = sourceRaster.getDataBufferIndex(x, y);
                    targetData.setElemFloatAt(index, masterData.getElemFloatAt(index));
                }
            }
        } catch (Throwable e) {
            OperatorUtils.catchOperatorException(getId(), e);
        } finally {
            progressMonitor.done();
        }
    }

    private synchronized void compute(final Band sourceBand,
            final ProgressMonitor progressMonitor) {

        if (processed) {
            return;
        }

        try {
            final RenderedImage renderedImage = sourceBand.getSourceImage().getImage(0);

            BufferedImage bufferedImage = new BufferedImage(sourceBand.getSceneRasterWidth(),
                    sourceBand.getSceneRasterHeight(),
                    BufferedImage.TYPE_USHORT_GRAY);
            bufferedImage.setData(PlanarImage.wrapRenderedImage(renderedImage).getAsBufferedImage().getRaster());

            final ImagePlus imagePlus = new ImagePlus(sourceBand.getDisplayName(), bufferedImage);
            ShortProcessor imageProcessor = new ShortProcessor(bufferedImage);
            imagePlus.setProcessor(imageProcessor);
            imagePlus.setCalibration(imagePlus.getCalibration());

            int width = imageProcessor.getWidth();
            int height = imageProcessor.getHeight();

            int intMax = (int) imageProcessor.getMax();
            System.out.println(intMax);

            N = width * height;
            probabilityHistogramDone = false;

            progressMonitor.worked(1);

            GrayLevel grayLevelTrue =
                    new GrayLevel(imageProcessor, true);
            GrayLevel grayLevelFalse =
                    new GrayLevel(imageProcessor, false);

            float fullMu = (float) grayLevelTrue.getOmega() * grayLevelTrue.getMu()
                    + (float) grayLevelFalse.getOmega() * grayLevelFalse.getMu();

            double sigmaMax = 0d;
            float threshold = 0f;

            for (int i = 0; i < intMax; i++) {

                double sigma = (double) grayLevelTrue.getOmega() * (Math.pow(grayLevelTrue.getMu() - fullMu, 2))
                        + (double) grayLevelFalse.getOmega()
                        * (Math.pow(grayLevelFalse.getMu() - fullMu, 2));

                if (sigma > sigmaMax) {
                    sigmaMax = sigma;
                    threshold = grayLevelTrue.getThreshold();
                }

                grayLevelTrue.addToEnd();
                grayLevelFalse.removeFromBeginning();
            }

            int offset = 0;

//            threshold = calculateThreshold(imageProcessor);

            short[] pixels = (short[]) imageProcessor.getPixels();

            for (int y = 0; y < height; y++) {
                offset = y * width;
                for (int x = 0; x < width; x++) {
                    int v = pixels[offset + x];
                    if (v <= threshold) {
                        imageProcessor.putPixel(x, y, 0);
                    } else {
                        imageProcessor.putPixel(x, y, 1);
                    }
                }
            }
//            File buffOutputFile = new File("D:/SOCIS/" + sourceBand.getDisplayName() + "-OtsuThreshold" + ".png");
//            try {
//                ImageIO.write(imageProcessor.get16BitBufferedImage(), "png", buffOutputFile);
//            } catch (IOException ex) {
//                Logger.getLogger(OtsuThresholdingOp.class.getName()).log(Level.SEVERE, null, ex);
//            }
//            IJ.setThreshold(threshold, intMax);
//                System.out.println();

//            imagePlus.setProcessor(imageProcessor.convertToShort(true));
//            imagePlus.setCalibration(imagePlus.getCalibration());

//            ImageConverter thresholdedImageConverter = new ImageConverter(imagePlus);
//            thresholdedImageConverter.convertToGray32();

            progressMonitor.worked(1);

//            RenderedImage thresholdedRenderedImage = toBufferedImage(imagePlus.getImage(),
//                    BufferedImage.TYPE_USHORT_GRAY);

            progressMonitor.worked(1);

            sourceBand.setSourceImage(imageProcessor.get16BitBufferedImage());

            progressMonitor.worked(1);

            processed = true;
        } catch (Throwable t) {
            OperatorUtils.catchOperatorException(this.getId(), t);
        } finally {
            progressMonitor.done();
        }
    }

    private static BufferedImage toBufferedImage(final Image image, final int type) {

        if (image instanceof BufferedImage) {
            return (BufferedImage) image;
        }
        if (image instanceof VolatileImage) {
            return ((VolatileImage) image).getSnapshot();
        }

        final BufferedImage buffImg = new BufferedImage(image.getWidth(null),
                image.getHeight(null), type);
        final Graphics2D g2 = buffImg.createGraphics();
//        g2.setRenderingHint(RenderingHints.KEY_INTERPOLATION,
//                RenderingHints.VALUE_INTERPOLATION_BICUBIC);
        g2.drawImage(image, null, null);
        g2.dispose();
        return buffImg;
    }
    private static final int NUM_BINS = 256;

    protected int[] makeHistogram(ImageProcessor imageProcessor) {
        int[] histData = new int[NUM_BINS];

        // Calculate histogram
        for (int r = 0; r < imageProcessor.getHeight(); r++) {
            for (int c = 0; c < imageProcessor.getWidth(); c++) {
                int h = (int) (imageProcessor.getPixel(r, c) * (NUM_BINS - 1));
                histData[h]++;
            }
        }

        return histData;
    }

    public float calculateThreshold(ImageProcessor imageProcessor) {

        int[] histData = imageProcessor.getHistogram();
        int total = imageProcessor.getWidth() * imageProcessor.getHeight();

        float sum = 0;
        for (int t = 0; t < NUM_BINS; t++) {
            sum += t * histData[t];
        }

        float sumB = 0;
        int wB = 0;
        int wF = 0;

        float varMax = 0;
        float threshold = 0f;

        for (int t = 0; t < NUM_BINS; t++) {
            wB += histData[t];               // Weight Background
            if (wB == 0) {
                continue;
            }

            wF = total - wB;                 // Weight Foreground
            if (wF == 0) {
                break;
            }

            sumB += (t * histData[t]);

            float mB = sumB / wB;            // Mean Background
            float mF = (sum - sumB) / wF;    // Mean Foreground

            // Calculate Between Class Variance
            float varBetween = (float) wB * (float) wF * (mB - mF) * (mB - mF);

            // Check if new maximum found
            if (varBetween > varMax) {
                varMax = varBetween;
                threshold = t;
            }
        }

        return (float) threshold / (NUM_BINS - 1);


    }

    private class GrayLevel {

        private int index;
        private float omega;
        private float mu;

        public GrayLevel(ImageProcessor imageProcessor, boolean isFirst) {

            if (!probabilityHistogramDone) {
                int[] histogram = imageProcessor.getHistogram();
                int length = histogram.length;
                probabilityHistogram = new float[length];

                for (int i = 0; i < length; i++) {
                    probabilityHistogram[i] = ((float) histogram[i]) / ((float) N);
                }
                probabilityHistogramDone = true;
            }

            if (isFirst) {
                index = 1;
                omega = probabilityHistogram[index - 1];
                if (omega == 0) {
                    mu = 0;
                } else {
                    mu = 1 * probabilityHistogram[index - 1] / omega;
                }
            } else {
                index = 2;
                omega = 0;
                mu = 0;
                for (int i = index; i < probabilityHistogram.length; i++) {
                    omega += probabilityHistogram[i - 1];
                    mu += probabilityHistogram[i - 1] * i;
                }
                if (omega == 0) {
                    mu = 0;
                } else {
                    mu /= omega;
                }
            }
        }

        public void removeFromBeginning() {
            index++;
            mu = 0;
            omega = 0;

            for (int i = index; i < probabilityHistogram.length; i++) {
                omega += probabilityHistogram[i - 1];
                mu += i * probabilityHistogram[i - 1];//i*
            }
            if (omega == 0) {
                mu = 0;
            } else {
                mu /= omega;
            }
        }

        public void addToEnd() {
            index++;
            mu = 0;
            omega = 0;
            for (int i = 1; i < index; i++) {
                omega += probabilityHistogram[i - 1];
                mu += i * probabilityHistogram[i - 1];
            }
            if (omega == 0) {
                mu = 0;
            } else {
                mu /= omega;
            }
        }

        public float getMu() {
            return mu;
        }

        public float getOmega() {
            return omega;
        }

        public int getThreshold() {
            return index;
        }
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
            super(OtsuThresholdingOp.class);
            setOperatorUI(OtsuThresholdingOpUI.class);
        }
    }
}
