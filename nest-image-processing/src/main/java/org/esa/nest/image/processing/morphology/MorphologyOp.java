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
package org.esa.nest.image.processing.morphology;

import com.bc.ceres.core.ProgressMonitor;
import ij.ImagePlus;
import ij.process.ImageConverter;
import ij.process.ImageProcessor;
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

import java.awt.*;
import java.awt.image.BufferedImage;
import java.awt.image.RenderedImage;
import java.awt.image.VolatileImage;
import java.util.HashMap;
import java.util.Map;
import javax.media.jai.PlanarImage;

/**
 * Morphology Operators (Dilate, Erode, Open, Close)
 */
@OperatorMetadata(alias = "MorphologyOperator",
category = "SAR Tools\\Image Processing",
description = "Morphology Operators")
public class MorphologyOp extends Operator {

    @SourceProduct(alias = "source")
    private Product sourceProduct = null;
    @TargetProduct
    private Product targetProduct;
    @Parameter(description = "The list of source bands.", alias = "sourceBands", itemAlias = "band",
    rasterDataNodeType = Band.class, label = "Source Bands")
    private String[] sourceBandNames;
    @Parameter(valueSet = {DILATE_OPERATOR, ERODE_OPERATOR, OPEN_OPERATOR,
        CLOSE_OPERATOR}, defaultValue = DILATE_OPERATOR,
    label = "Operator")
    private String operator;
    @Parameter(description = "Iterations", interval = "[1, 100]", defaultValue = "1", label = "Iterations")
    private int nIterations = 1;
    static final String DILATE_OPERATOR = "Dilate";
    static final String ERODE_OPERATOR = "Erode";
    static final String OPEN_OPERATOR = "Open";
    static final String CLOSE_OPERATOR = "Close";
    private final Map<String, String[]> targetBandNameToSourceBandName = new HashMap<String, String[]>();
    private int sourceImageWidth;
    private int sourceImageHeight;
    private boolean processed = false;

    public MorphologyOp() {
    }

    /**
     * Set speckle filter. This function is used by unit test only.
     *
     * @param s The filter name.
     */
    public void SetFilter(String s) {

        if (s.equals(DILATE_OPERATOR)
                || s.equals(ERODE_OPERATOR)
                || s.equals(OPEN_OPERATOR)
                || s.equals(CLOSE_OPERATOR)) {
            operator = s;
        } else {
            throw new OperatorException(s + " is an invalid filter name.");
        }
    }

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

            Tile sourceRaster;
            final String[] srcBandNames = targetBandNameToSourceBandName.get(targetBand.getName());
            Band sourceBand = sourceProduct.getBand(srcBandNames[0]);
            sourceRaster = getSourceTile(sourceBand, targetTileRectangle);

            if (sourceRaster == null) {
                throw new OperatorException("Cannot get source tile");
            }
            if (!processed) {
                compute(sourceBand, pm);
            }
            final ProductData targetData = targetTile.getDataBuffer();
            final ProductData sourceData = sourceRaster.getDataBuffer();

            for (int y = y0; y < y0 + h; y++) {
                for (int x = x0; x < x0 + w; x++) {
                    final int index = sourceRaster.getDataBufferIndex(x, y);
                    targetData.setElemFloatAt(targetTile.getDataBufferIndex(x, y), sourceData.getElemFloatAt(index));
                }
            }
        } catch (Throwable e) {
            OperatorUtils.catchOperatorException(getId(), e);
        } finally {
            pm.done();
        }
    }

    private synchronized void compute(final Band sourceBand,
            final ProgressMonitor progressMonitor) {

        if (processed) {
            return;
        }

        try {
            final RenderedImage ri = sourceBand.getSourceImage().getImage(0);

            BufferedImage img = new BufferedImage(sourceBand.getSceneRasterWidth(), sourceBand.getSceneRasterHeight(),
                    BufferedImage.TYPE_USHORT_GRAY);
            img = PlanarImage.wrapRenderedImage(ri).getAsBufferedImage();
            final ImagePlus imp = new ImagePlus(sourceBand.getDisplayName(), img);
            ImageProcessor ip = imp.getProcessor();

            for (int i = 0; i < nIterations; i++) {
                if (operator.equals(DILATE_OPERATOR)) {
                    ip.dilate();
                } else if (operator.equals(ERODE_OPERATOR)) {
                    ip.erode();
                } else if (operator.equals(CLOSE_OPERATOR)) {
                    ip.dilate();
                    progressMonitor.worked(1);
                    ip.erode();
                } else if (operator.equals(OPEN_OPERATOR)) {
                    ip.erode();
                    progressMonitor.worked(1);
                    ip.dilate();
                }
            }

            progressMonitor.worked(1);

            imp.setProcessor(ip.convertToShort(true));
            imp.setCalibration(imp.getCalibration());
            ImageConverter ic = new ImageConverter(imp);
            ic.convertToGray32();

            progressMonitor.worked(1);

            RenderedImage renderedImage = toBufferedImage(imp.getImage(), BufferedImage.TYPE_USHORT_GRAY);

            progressMonitor.worked(1);

            sourceBand.setSourceImage(renderedImage);

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
        g2.setRenderingHint(RenderingHints.KEY_INTERPOLATION,
                RenderingHints.VALUE_INTERPOLATION_BICUBIC);
        g2.drawImage(image, null, null);
        g2.dispose();
        return buffImg;
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
            super(MorphologyOp.class);
            setOperatorUI(MorphologyOpUI.class);
        }
    }
}
