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
import java.net.URL;
import java.util.Arrays;
import org.esa.beam.dataio.envisat.EnvisatProductReader;
import org.esa.beam.dataio.envisat.EnvisatProductReaderPlugIn;
import org.esa.beam.framework.datamodel.Band;
import org.esa.beam.framework.datamodel.Product;
import org.esa.beam.framework.gpf.OperatorSpi;
import org.esa.nest.image.processing.segmentation.thresholding.BasicThresholdingOp;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author Emanuela Boros
 * @since October 2012
 */
public class MaximumEntropyThresholdingOpTest {

    private OperatorSpi spi;
    private final static String inputSARpath = "SAR/"
            + "ASA_WSM_1PNUPA20080124_101559_000000672065_00237_30852_1983.N1";
    EnvisatProductReaderPlugIn readerPlugIn = new EnvisatProductReaderPlugIn();
    EnvisatProductReader reader;
    Product sourceProduct;
    MaximumEntropyThresholdingOp op;
    

    @Before
    public void setUp() throws Exception {
        URL url = getClass().getClassLoader().getResource(inputSARpath);
        spi = new MaximumEntropyThresholdingOp.Spi();
        reader = (EnvisatProductReader) readerPlugIn.createReaderInstance();
        sourceProduct = reader.readProductNodes(url.getFile(), null);
        op = (MaximumEntropyThresholdingOp) spi.createOperator();
    }

    @After
    public void tearDown() {
        spi = null;
        readerPlugIn = null;
    }

    /**
     * Tests MaximumEntropyThresholding with a 4-by-4 test product.
     *
     * @throws Exception The exception.
     */
    @Test
    public void testMaximumEntropyThresholdingOperator() throws Exception {

        assertNotNull(sourceProduct);
        assertNotNull(op);
        op.setSourceProduct(sourceProduct);

        final Product targetProduct = op.getTargetProduct();
        final Band band = targetProduct.getBandAt(0);
        assertNotNull(band);

        final float[] floatValues = new float[8];
        band.readPixels(345, 300, 4, 2, floatValues, ProgressMonitor.NULL);
        for (float f : floatValues) {
            System.out.print(f + ",");
        }
        final float[] expected = {-1f, 0f, -1f, -1f, -1f, -1f, -1f, -1f};
        assertTrue(Arrays.equals(expected, floatValues));
    }
}