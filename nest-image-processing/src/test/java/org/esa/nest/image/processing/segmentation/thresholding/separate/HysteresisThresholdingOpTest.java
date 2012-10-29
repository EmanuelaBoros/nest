/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
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
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author Emanuela
 */
public class HysteresisThresholdingOpTest {

    private OperatorSpi spi;
    private final static String inputSARpath = "SAR/"
            + "ASA_WSM_1PNUPA20080124_101559_000000672065_00237_30852_1983.N1";
    EnvisatProductReaderPlugIn readerPlugIn = new EnvisatProductReaderPlugIn();
    EnvisatProductReader reader;
    Product sourceProduct;
    HysteresisThresholdingOp op;
    

    @Before
    public void setUp() throws Exception {
        URL url = getClass().getClassLoader().getResource(inputSARpath);
        spi = new HysteresisThresholdingOp.Spi();
        reader = (EnvisatProductReader) readerPlugIn.createReaderInstance();
        sourceProduct = reader.readProductNodes(url.getFile(), null);
        op = (HysteresisThresholdingOp) spi.createOperator();
    }

    @After
    public void tearDown() {
        spi = null;
        readerPlugIn = null;
    }

    /**
     * Tests HysteresisThresholdingOperator with a 4-by-4 test product.
     *
     * @throws Exception The exception.
     */
    @Test
    public void testHysteresisThresholdingOperator() throws Exception {

        assertNotNull(sourceProduct);
        assertNotNull(op);
        op.setSourceProduct(sourceProduct);

        final Product targetProduct = op.getTargetProduct();
        final Band band = targetProduct.getBandAt(0);
        assertNotNull(band);

        final float[] floatValues = new float[8];
        band.readPixels(58, 184, 4, 2, floatValues, ProgressMonitor.NULL);

        final float[] expected = {0f, -1f, -1f, -1f, 0f, -1f, -1f, -1f};
        assertTrue(Arrays.equals(expected, floatValues));
    }
}
