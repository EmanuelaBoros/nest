package org.esa.nest.gpf;

import junit.framework.TestCase;
import org.esa.beam.framework.gpf.OperatorSpi;
import org.esa.beam.framework.gpf.GPF;
import org.esa.beam.framework.datamodel.*;
import org.esa.nest.datamodel.AbstractMetadata;
import org.esa.nest.gpf.TestOperator;

import java.util.Arrays;
import java.io.File;

import com.bc.ceres.core.ProgressMonitor;

/**
 * Unit test for UndersamplingOperator.
 */
public class TestUndersamplingOperator extends TestCase {

    private OperatorSpi spi;

    @Override
    protected void setUp() throws Exception {
        spi = new UndersamplingOp.Spi();
        GPF.getDefaultInstance().getOperatorSpiRegistry().addOperatorSpi(spi);
    }

    @Override
    protected void tearDown() throws Exception {
        GPF.getDefaultInstance().getOperatorSpiRegistry().removeOperatorSpi(spi);
    }

    /**
     * Tests undersampling operator with a 6x12 "DETECTED" test product.
     * @throws Exception general exception
     */
    public void testUndersampling() throws Exception {

        Product sourceProduct = createTestProduct(12, 6);

        UndersamplingOp op = (UndersamplingOp)spi.createOperator();
        assertNotNull(op);
        op.setSourceProduct(sourceProduct);

        op.setUndersamplingMethod("Kernel Filtering");
        op.setFilterType("Loss Pass");
        op.setOutputImageSaize(2, 4);

        String home = "C:" + File.separator + "NEST" + File.separator + "nest" + File.separator + "beam";
        System.setProperty("nest.home", home);

        // get targetProduct: execute initialize()
        Product targetProduct = op.getTargetProduct();
        TestOperator.verifyProduct(targetProduct);

        Band band = targetProduct.getBandAt(0);
        assertNotNull(band);

        // readPixels: execute computeTiles()
        float[] floatValues = new float[8];
        band.readPixels(0, 0, 4, 2, floatValues, ProgressMonitor.NULL);

        // compare with expected outputs:
        float[] expectedValues = {126.0f, 153.0f,   180.0f, 207.0f, 450.0f, 477.0f, 504.0f, 531.0f};
        assertTrue(Arrays.equals(expectedValues, floatValues));

        // compare updated metadata

        MetadataElement abs = targetProduct.getMetadataRoot().getElement(Product.ABSTRACTED_METADATA_ROOT_NAME);

        TestOperator.attributeEquals(abs, AbstractMetadata.azimuth_spacing, 4.5);
        TestOperator.attributeEquals(abs, AbstractMetadata.range_spacing, 6.0);
        TestOperator.attributeEquals(abs, AbstractMetadata.line_time_interval, 0.03);
        TestOperator.attributeEquals(abs, AbstractMetadata.first_line_time, "10-MAY-2008 20:32:46.895683");        
    }


    /**
     * Creates a 6-by-12 test product as shown below:
     *  1  2  3  4  5  6  7  8  9 10 11 12
     * 13 14 15 16 17 18 19 20 21 22 23 24
     * 25 26 27 28 29 30 31 32 33 34 35 36
     * 37 38 39 40 41 42 43 44 45 46 47 48
     * 49 50 51 52 53 54 55 56 57 58 59 60
     * 61 62 63 64 65 66 67 68 69 70 71 72
     * @param w width
     * @param h height
     * @return the created product
     */
    private static Product createTestProduct(int w, int h) {

        Product testProduct = TestOperator.createProduct("ASA_APG_1P", w, h);

        // create a Band: band1
        Band band1 = testProduct.addBand("band1", ProductData.TYPE_INT32);
        band1.setUnit("amplitude");
        band1.setSynthetic(true);
        int[] intValues = new int[w * h];
        for (int i = 0; i < w * h; i++) {
            intValues[i] = i + 1;
        }
        band1.setData(ProductData.createInstance(intValues));

        // create abstracted metadata
        MetadataElement abs = testProduct.getMetadataRoot().getElement(Product.ABSTRACTED_METADATA_ROOT_NAME);

        AbstractMetadata.setAttribute(abs, AbstractMetadata.range_spacing, 2.0F);
        AbstractMetadata.setAttribute(abs, AbstractMetadata.azimuth_spacing, 1.5F);
        AbstractMetadata.setAttribute(abs, AbstractMetadata.line_time_interval, 0.01F);
        AbstractMetadata.setAttribute(abs, AbstractMetadata.first_line_time,
                AbstractMetadata.parseUTC("10-MAY-2008 20:32:46.885684"));

        return testProduct;
    }
}