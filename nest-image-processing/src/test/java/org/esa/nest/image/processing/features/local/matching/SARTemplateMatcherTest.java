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
package org.esa.nest.image.processing.features.local.matching;

import java.awt.Rectangle;
import java.net.URL;
import org.esa.beam.dataio.envisat.EnvisatProductReader;
import org.esa.beam.dataio.envisat.EnvisatProductReaderPlugIn;
import org.esa.beam.framework.datamodel.Band;
import org.esa.beam.framework.datamodel.Product;
import org.esa.beam.framework.gpf.OperatorSpi;
import org.esa.beam.framework.gpf.Tile;
import org.esa.nest.image.processing.features.local.matching.SARTemplateMatcher.Mode;
import org.esa.nest.image.processing.segmentation.thresholding.BasicThresholdingOp;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author Emanuela
 */
public class SARTemplateMatcherTest {

    private OperatorSpi spi;
    private final static String inputSARpath = "SAR/"
            + "ASA_WSM_1PNUPA20080124_101559_000000672065_00237_30852_1983.N1";
    EnvisatProductReaderPlugIn readerPlugIn = new EnvisatProductReaderPlugIn();
    private int sourceImageWidth;
    private int sourceImageHeight;
    private int halfSizeX;
    private int halfSizeY;
    private int filterSizeX = 3;
    private int filterSizeY = 3;

    @Before
    public void setUp() {
        spi = new BasicThresholdingOp.Spi();
    }

    @After
    public void tearDown() {
        spi = null;
        readerPlugIn = null;
    }

    /**
     * Test of computeMatchScore method, of class SARTemplateMatcher.
     */
    @Test
    public void testComputeMatchScore() throws Exception {

        final EnvisatProductReader reader = (EnvisatProductReader) readerPlugIn.createReaderInstance();
        URL url = getClass().getClassLoader().getResource(inputSARpath);
        final Product sourceProduct = reader.readProductNodes(url.getFile(), null);
        assertNotNull(sourceProduct);
        final BasicThresholdingOp op = (BasicThresholdingOp) spi.createOperator();
        assertNotNull(op);
        op.setSourceProduct(sourceProduct);
        Band sourceBand = sourceProduct.getBandAt(0);

        sourceImageWidth = sourceProduct.getSceneRasterWidth();
        sourceImageHeight = sourceProduct.getSceneRasterHeight();

        halfSizeX = filterSizeX / 2;
        halfSizeY = filterSizeY / 2;

        final Rectangle sourceTileRectangle = getSourceTileRectangle(58, 184, 4, 2);
        final Rectangle templateRectangle = getSourceTileRectangle(58, 184, 4, 2);

        Tile sourceTile = op.getSourceTile(sourceBand, sourceTileRectangle);
        Tile templateTile = op.getSourceTile(sourceBand, templateRectangle);

        SARTemplateMatcher instance = new SARTemplateMatcher();
        float expResult = 0.0F;
        float result = instance.computeMatchScore(sourceBand, sourceTile, templateTile,
                Mode.SUM_SQUARED_DIFFERENCE);
        assertEquals(expResult, result, 0f);
    }

    /**
     * Test of computeMatchScore method, of class SIFTMatcher.
     */
    @Test
    public void testComputeSIFTMatchScore() throws Exception {

        final EnvisatProductReader reader = (EnvisatProductReader) readerPlugIn.createReaderInstance();
        URL url = getClass().getClassLoader().getResource(inputSARpath);
        final Product sourceProduct = reader.readProductNodes(url.getFile(), null);
        assertNotNull(sourceProduct);
        final BasicThresholdingOp op = (BasicThresholdingOp) spi.createOperator();
        assertNotNull(op);
        op.setSourceProduct(sourceProduct);
        Band sourceBand = sourceProduct.getBandAt(0);

        sourceImageWidth = sourceProduct.getSceneRasterWidth();
        sourceImageHeight = sourceProduct.getSceneRasterHeight();

        halfSizeX = filterSizeX / 2;
        halfSizeY = filterSizeY / 2;

        final Rectangle sourceTileRectangle = getSourceTileRectangle(58, 184, 4, 2);
        final Rectangle templateRectangle = getSourceTileRectangle(100, 245, 4, 2);

        Tile sourceTile = op.getSourceTile(sourceBand, sourceTileRectangle);
        Tile templateTile = op.getSourceTile(sourceBand, templateRectangle);

        SIFTMatcher instance = new SIFTMatcher();
        float fExpResult = 0f;
        float fResult = instance.computeMatchScore(sourceBand, sourceTile, templateTile,
                Mode.SUM_SQUARED_DIFFERENCE);
        System.out.println("Result: " + fResult);
        assertEquals(fExpResult, fResult, 0f);
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
