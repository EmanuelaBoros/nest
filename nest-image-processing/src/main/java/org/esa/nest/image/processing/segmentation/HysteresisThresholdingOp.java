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
import ij.ImageStack;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;
import java.awt.Rectangle;
import java.awt.image.BufferedImage;
import java.awt.image.RenderedImage;
import java.util.HashMap;
import java.util.Map;
import org.esa.beam.framework.datamodel.Band;
import org.esa.beam.framework.datamodel.Product;
import org.esa.beam.framework.datamodel.ProductData;
import org.esa.beam.framework.gpf.Operator;
import org.esa.beam.framework.gpf.OperatorException;
import org.esa.beam.framework.gpf.Tile;
import org.esa.beam.framework.gpf.annotations.OperatorMetadata;
import org.esa.beam.framework.gpf.annotations.Parameter;
import org.esa.beam.framework.gpf.annotations.SourceProduct;
import org.esa.beam.framework.gpf.annotations.TargetProduct;
import org.esa.nest.gpf.OperatorUtils;

@OperatorMetadata(alias = "HysteresisThresholding",
category = "SAR Tools\\Image Processing",
description = "HysteresisThresholding")
public class HysteresisThresholdingOp extends Operator {

    public static float[] probabilityHistogram;
    private float threshold;
    final static int MAX_VALUE = 0;
    final static int MIN_VALUE = 256;
    private boolean probabilityHistogramDone;
    public static int N;
    @SourceProduct(alias = "source")
    private Product sourceProduct = null;
    @TargetProduct
    private Product targetProduct;
    @Parameter(description = "The list of source bands.", alias = "sourceBands", itemAlias = "band",
    rasterDataNodeType = Band.class, label = "Source Bands")
    private String[] sourceBandNames;
    @Parameter(description = "HighThreshold", defaultValue = "5.0", label = "HighThreshold")
    private double highThreshold = 5.0;
    @Parameter(description = "LowThreshold", defaultValue = "1.0", label = "LowThreshold")
    private double lowThreshold = 1.0;
    private final Map<String, String[]> targetBandNameToSourceBandName = new HashMap<String, String[]>();
    private int sourceImageWidth;
    private int sourceImageHeight;
    private boolean processed = false;
    private int halfSizeX;
    private int halfSizeY;
    private int filterSizeX = 3;
    private int filterSizeY = 3;
    private static ImagePlus fullImagePlus;

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
     * during operator initialisation.
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

            computeOtsuThesholding(sourceBand, sourceRaster, targetTile, x0, y0, w, h, pm);

        } catch (Throwable e) {
            OperatorUtils.catchOperatorException(getId(), e);
        } finally {
            pm.done();
        }
    }

    /**
     * Apply Otsu Thresholding
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
     * computation cancelation requests.
     * @throws org.esa.beam.framework.gpf.OperatorException If an error occurs
     * during computation of the filtered value.
     */
    private void computeOtsuThesholding(final Band sourceBand, final Tile sourceRaster,
            final Tile targetTile, final int x0, final int y0, final int w, final int h,
            final ProgressMonitor pm) {

        if (!processed) {
            final RenderedImage fullRenderedImage = sourceBand.getSourceImage().getImage(0);
            final BufferedImage fullBufferedImage = new BufferedImage(sourceBand.getSceneRasterWidth(),
                    sourceBand.getSceneRasterHeight(),
                    BufferedImage.TYPE_USHORT_GRAY);
            fullBufferedImage.setData(fullRenderedImage.getData());

            fullImagePlus = new ImagePlus(sourceBand.getDisplayName(), fullBufferedImage);
            processed = true;
        }

        final ImageProcessor fullImageProcessor = fullImagePlus.getProcessor();

        ByteProcessor fullByteProcessor = (ByteProcessor) fullImageProcessor.convertToByte(true);

        int width = fullByteProcessor.getWidth();
        int height = fullByteProcessor.getHeight();

        int intMax = (int) fullByteProcessor.getMax();
        N = width * height;
        probabilityHistogramDone = false;

        ImageStack stack = fullImagePlus.getStack();
        ImageStack res_trin = new ImageStack(stack.getWidth(), stack.getHeight());
        ImageStack res_hyst = new ImageStack(stack.getWidth(), stack.getHeight());
        ImageProcessor tmp1;
        ImageProcessor tmp2;

        for (int s = 1; s <= stack.getSize(); s++) {
            tmp1 = trinarise(stack.getProcessor(s), highThreshold, lowThreshold);
            tmp2 = hysteresisThresholding(tmp1);
            res_trin.addSlice("", tmp1);
            res_hyst.addSlice("", tmp2);
        }

        final Rectangle srcTileRectangle = sourceRaster.getRectangle();

        fullByteProcessor.setRoi(srcTileRectangle);

        ImageProcessor roiImageProcessor = fullByteProcessor.crop();

        int offset = 0;

        byte[] pixels = (byte[]) roiImageProcessor.getPixels();

        for (int y = 0; y < roiImageProcessor.getHeight(); y++) {
            offset = y * roiImageProcessor.getWidth();
            for (int x = 0; x < roiImageProcessor.getWidth(); x++) {
                int value = pixels[offset + x];
                if (value > threshold) {
                    roiImageProcessor.putPixel(x, y, intMax);
                } else {
                    roiImageProcessor.putPixel(x, y, 0);
                }
            }
        }

        final ProductData trgData = targetTile.getDataBuffer();
        final ProductData sourceData = ProductData.createInstance((byte[]) roiImageProcessor.getPixels());

        final int maxY = y0 + h;
        final int maxX = x0 + w;
        for (int y = y0; y < maxY; ++y) {
            for (int x = x0; x < maxX; ++x) {

                trgData.setElemDoubleAt(targetTile.getDataBufferIndex(x, y),
                        sourceData.getElemDoubleAt(sourceRaster.getDataBufferIndex(x, y)));
            }
        }
    }

    /**
     * Double thresholding
     *
     * @param imageProcessor original image
     * @param highThreshold high threshold
     * @param lowThreshold low threshold
     * @return "trinarised" image
     */
    ImageProcessor trinarise(ImageProcessor imageProcessor, double highThreshold,
            double lowThreshold) {

        int la = imageProcessor.getWidth();
        int ha = imageProcessor.getHeight();
        ByteProcessor res = new ByteProcessor(la, ha);
        float pix;

        for (int x = 0; x < la; x++) {
            for (int y = 0; y < ha; y++) {
                pix = imageProcessor.getPixelValue(x, y);
                if (pix >= highThreshold) {
                    res.putPixel(x, y, 255);
                } else if (pix >= lowThreshold) {
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
    ImageProcessor hysteresisThresholding(ImageProcessor ima) {
        int la = ima.getWidth();
        int ha = ima.getHeight();
        ImageProcessor res = ima.duplicate();
        float pix;
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
}
