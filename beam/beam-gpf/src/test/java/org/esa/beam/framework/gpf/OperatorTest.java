package org.esa.beam.framework.gpf;

import com.bc.ceres.core.ProgressMonitor;
import junit.framework.TestCase;
import org.esa.beam.framework.datamodel.Band;
import org.esa.beam.framework.datamodel.Product;
import org.esa.beam.framework.datamodel.ProductData;
import org.esa.beam.framework.gpf.annotations.TargetProduct;
import org.esa.beam.framework.gpf.operators.common.PassThroughOp;

import java.io.IOException;


public class OperatorTest extends TestCase {


    public void testBasicOperatorStates() throws OperatorException, IOException {
        final FooOp op = new FooOp();
        assertNotNull(op.getSpi());
        assertFalse(op.initializeCalled);
        assertFalse(op.computeTileCalled);
        final Product product = op.getTargetProduct();
        assertNotNull(product);
        assertFalse("products created by operators cannot be modified", product.isModified());
        assertTrue("op.initialize not called", op.initializeCalled);
        assertFalse("op.computeTileCalled called", op.computeTileCalled);
        product.getBand("bar").readRasterDataFully(ProgressMonitor.NULL);
        assertTrue(op.initializeCalled);
        assertTrue(op.computeTileCalled);
    }

    public void testPassThroughDetection() throws OperatorException, IOException {
        Product sourceProduct = createFooProduct();
        final Operator op = new PassThroughOp(sourceProduct);
        assertNotNull(op.getSpi());
        assertFalse(op.context.isPassThrough());
        Product targetProduct = op.getTargetProduct();// force init
        assertSame(sourceProduct, targetProduct);
        assertTrue(op.context.isPassThrough());
    }

    public void testSourceProducts() throws IOException, OperatorException {
        final Operator operator = new Operator() {
            @Override
            public void initialize() throws OperatorException {
            }
        };

        final Product sp1 = new Product("sp1", "t", 1, 1);
        final Product sp2 = new Product("sp2", "t", 1, 1);
        final Product sp3 = new Product("sp3", "t", 1, 1);

        operator.setSourceProduct(sp1);
        assertSame(sp1, operator.getSourceProduct());
        assertSame(sp1, operator.getSourceProduct("sourceProduct"));

        operator.setSourceProduct("sp1", sp1);
        assertSame(sp1, operator.getSourceProduct());
        assertSame(sp1, operator.getSourceProduct("sourceProduct"));
        assertSame(sp1, operator.getSourceProduct("sp1"));

        Product[] products = operator.getSourceProducts();
        assertNotNull(products);
        assertEquals(1, products.length);
        assertSame(sp1, products[0]);

        operator.setSourceProduct("sp2", sp2);
        products = operator.getSourceProducts();
        assertNotNull(products);
        assertEquals(2, products.length);
        assertSame(sp1, products[0]);
        assertSame(sp2, products[1]);

        operator.setSourceProducts(new Product[]{sp3, sp2, sp1});
        assertNull(operator.getSourceProduct("sourceProduct"));
        assertNull(operator.getSourceProduct("sp1"));
        assertNull(operator.getSourceProduct("sp2"));
        products = operator.getSourceProducts();
        assertNotNull(products);
        assertEquals(3, products.length);
        assertSame(sp3, products[0]);
        assertSame(sp2, products[1]);
        assertSame(sp1, products[2]);
        assertSame(sp3, operator.getSourceProduct("sourceProduct1"));
        assertSame(sp2, operator.getSourceProduct("sourceProduct2"));
        assertSame(sp1, operator.getSourceProduct("sourceProduct3"));
        assertEquals("sourceProduct3", operator.getSourceProductId(sp1));
        assertEquals("sourceProduct2", operator.getSourceProductId(sp2));
        assertEquals("sourceProduct1", operator.getSourceProductId(sp3));


        operator.setSourceProducts(new Product[]{sp1, sp2, sp1});
        products = operator.getSourceProducts();
        assertNotNull(products);
        assertEquals(2, products.length);
        assertSame(sp1, products[0]);
        assertSame(sp2, products[1]);
        assertSame(sp1, operator.getSourceProduct("sourceProduct1"));
        assertSame(sp2, operator.getSourceProduct("sourceProduct2"));
        assertSame(sp1, operator.getSourceProduct("sourceProduct3"));
        assertEquals("sourceProduct1", operator.getSourceProductId(sp1));
        assertEquals("sourceProduct2", operator.getSourceProductId(sp2));
        assertNull(operator.getSourceProductId(sp3));
    }

    private static Product createFooProduct() {
        Product product = new Product("foo", "grunt", 1, 1);
        product.addBand("bar", ProductData.TYPE_FLOAT64);
        return product;
    }

    private static class FooOp extends Operator {
        private boolean initializeCalled;
        private boolean computeTileCalled;
        @TargetProduct
        private Product targetProduct;

        @Override
        public void initialize() throws OperatorException {
            initializeCalled = true;
            targetProduct = createFooProduct();
            targetProduct.addBand("foo", ProductData.TYPE_INT8); // will set the "modified" flag
        }

        @Override
        public void computeTile(Band band, Tile tile, ProgressMonitor pm) throws OperatorException {
            computeTileCalled = true;
        }
    }
}
