package org.esa.beam.framework.gpf;

import com.bc.ceres.core.ProgressMonitor;
import org.esa.beam.framework.datamodel.Band;
import org.esa.beam.framework.datamodel.Product;
import org.esa.beam.framework.datamodel.ProductData;
import org.esa.beam.framework.gpf.annotations.*;

import java.awt.Rectangle;
import java.util.Map;

public class TestOps {

    public static final int RASTER_WIDTH = 3;
    public static final int RASTER_HEIGHT = 2;
    static String calls = "";

    public static void registerCall(String str) {
        calls += str;
    }

    public static void clearCalls() {
        calls = "";
    }

    public static String getCalls() {
        return calls;
    }

    @OperatorMetadata(alias = "Op1")
    public static class Op1 extends Operator {
        @TargetProduct
        private Product targetProduct;

        @Override
        public void initialize() {
            targetProduct = new Product("Op1Name", "Op1Type", RASTER_WIDTH, RASTER_HEIGHT);
            targetProduct.addBand(new Band("Op1A", ProductData.TYPE_INT8, RASTER_WIDTH, RASTER_HEIGHT));
        }

        @Override
        public void computeTile(Band band, Tile targetTile, ProgressMonitor pm) {
            //System.out.println("=====>>>>>> Op1.computeBand  start");
            registerCall("Op1;");
            //System.out.println("=====>>>>>> Op1.computeBand  end");
        }

        public static class Spi extends OperatorSpi {

            public Spi() {
                super(Op1.class);
            }
        }
    }

    @OperatorMetadata(alias = "Op2")
    public static class Op2 extends Operator {

        @Parameter
        public double threshold;

        @SourceProduct(bands = {"Op1A"})
        public Product input;

        @TargetProduct
        public Product output;

        @Override
        public void initialize() {
            output = new Product("Op2Name", "Op2Type", RASTER_WIDTH, RASTER_HEIGHT);
            output.addBand(new Band("Op2A", ProductData.TYPE_INT8, RASTER_WIDTH, RASTER_HEIGHT));
            output.addBand(new Band("Op2B", ProductData.TYPE_INT8, RASTER_WIDTH, RASTER_HEIGHT));
        }

        @Override
        public void computeTileStack(Map<Band, Tile> targetTiles, Rectangle rectangle, ProgressMonitor pm) throws OperatorException {
            //System.out.println("=====>>>>>> Op2.computeAllBands  start");
            Tile tile1A = getSourceTile(input.getBand("Op1A"), rectangle, pm);

            Tile tile2A = targetTiles.get(output.getBand("Op2A"));
            Tile tile2B = targetTiles.get(output.getBand("Op2B"));
            //System.out.println("=====>>>>>> Op2.computeAllBands end");

            registerCall("Op2;");
        }

        public static class Spi extends OperatorSpi {
            public Spi() {
                super(Op2.class);
            }
        }
    }

    @OperatorMetadata(alias = "Op3")
    public static class Op3 extends Operator {

        @Parameter
        public boolean ignoreSign;

        @Parameter(description = "The valid mask expression")
        public String expression;

        @Parameter(valueSet = {"NN", "BQ", "CC"}, defaultValue = "NN")
        public String interpolMethod;

        @Parameter(defaultValue = "1.5", interval = "[-10,10)")
        public double factor;

        @SourceProduct(bands = {"Op1A"})
        public Product input1;

        @SourceProduct(bands = {"Op2A", "Op2B"})
        public Product input2;

        @SourceProducts
        public Product[] inputs;

        @TargetProduct
        public Product output;

        @Override
        public void initialize() {
            output = new Product("Op3Name", "Op3Type", RASTER_WIDTH, RASTER_HEIGHT);
            output.addBand(new Band("Op3A", ProductData.TYPE_INT8, RASTER_WIDTH, RASTER_HEIGHT));
            output.addBand(new Band("Op3B", ProductData.TYPE_INT8, RASTER_WIDTH, RASTER_HEIGHT));
            output.addBand(new Band("Op3C", ProductData.TYPE_INT8, RASTER_WIDTH, RASTER_HEIGHT));
            output.addBand(new Band("Op3D", ProductData.TYPE_INT8, RASTER_WIDTH, RASTER_HEIGHT));
        }

        @Override
        public void computeTileStack(Map<Band, Tile> targetTiles, Rectangle rectangle, ProgressMonitor pm) throws OperatorException {
            //System.out.println("=====>>>>>> Op3.computeAllBands  start");

            Tile tile1A = getSourceTile(input1.getBand("Op1A"), rectangle, pm);
            Tile tile2A = getSourceTile(input2.getBand("Op2A"), rectangle, pm);
            Tile tile2B = getSourceTile(input2.getBand("Op2B"), rectangle, pm);

            Tile tile3A = targetTiles.get(output.getBand("Op3A"));
            Tile tile3B = targetTiles.get(output.getBand("Op3B"));
            Tile tile3C = targetTiles.get(output.getBand("Op3C"));
            Tile tile3D = targetTiles.get(output.getBand("Op3D"));
            registerCall("Op3;");

            //System.out.println("=====>>>>>> Op3.computeAllBands  end");
        }

        public static class Spi extends OperatorSpi {
            public Spi() {
                super(Op3.class);
            }
        }
    }

    @OperatorMetadata(alias = "Op4")
    public static class Op4 extends Operator {
        @TargetProduct
        private Product targetProduct;

        @TargetProperty(alias = "PI", description = "The ratio of any circle's circumference to its diameter")
        private double pi;
        
        @TargetProperty
        private String[] names;

        @Override
        public void initialize() {
            targetProduct = new Product("Op1Name", "Op1Type", RASTER_WIDTH, RASTER_HEIGHT);
            targetProduct.addBand(new Band("Op1A", ProductData.TYPE_INT8, RASTER_WIDTH, RASTER_HEIGHT));
            pi = 3.142;
        }

        @Override
        public void computeTile(Band band, Tile targetTile, ProgressMonitor pm) {
            //System.out.println("=====>>>>>> Op4.computeBand  start");
            registerCall("Op4;");
            //System.out.println("=====>>>>>> Op4.computeBand  end");
        }

        public static class Spi extends OperatorSpi {

            public Spi() {
                super(Op4.class);
            }
        }
    }
}