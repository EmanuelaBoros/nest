/*
 * $Id: BandMathsOp.java,v 1.1 2010-03-31 13:59:24 lveci Exp $
 *
 * Copyright (C) 2007 by Brockmann Consult (info@brockmann-consult.de)
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation. This program is distributed in the hope it will
 * be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */
package org.esa.beam.gpf.operators.standard;

import com.bc.ceres.core.ProgressMonitor;
import com.bc.jexp.Namespace;
import com.bc.jexp.ParseException;
import com.bc.jexp.Parser;
import com.bc.jexp.Symbol;
import com.bc.jexp.Term;
import com.bc.jexp.WritableNamespace;
import com.bc.jexp.impl.ParserImpl;
import com.bc.jexp.impl.SymbolFactory;
import org.esa.beam.framework.datamodel.Band;
import org.esa.beam.framework.datamodel.Product;
import org.esa.beam.framework.datamodel.ProductData;
import org.esa.beam.framework.dataop.barithm.BandArithmetic;
import org.esa.beam.framework.dataop.barithm.BandArithmetic.ProductPrefixProvider;
import org.esa.beam.framework.dataop.barithm.RasterDataEvalEnv;
import org.esa.beam.framework.dataop.barithm.RasterDataSymbol;
import org.esa.beam.framework.gpf.Operator;
import org.esa.beam.framework.gpf.OperatorException;
import org.esa.beam.framework.gpf.OperatorSpi;
import org.esa.beam.framework.gpf.Tile;
import org.esa.beam.framework.gpf.ui.BandArithmeticOpUI;
import org.esa.beam.framework.gpf.annotations.OperatorMetadata;
import org.esa.beam.framework.gpf.annotations.Parameter;
import org.esa.beam.framework.gpf.annotations.SourceProducts;
import org.esa.beam.framework.gpf.annotations.TargetProduct;
import org.esa.beam.util.ProductUtils;
import org.esa.beam.util.StringUtils;

import java.awt.Rectangle;
import java.util.HashMap;
import java.util.Map;

@OperatorMetadata(alias = "BandMath", version = "1.0", category = "Utilities",
                  description = "Create a product with one or more bands using mathematical expressions.\n" +
                                "This operator can only be invoked with a Graph XML file.", internal=true)
public class BandMathsOp extends Operator {

    public static class BandDescriptor {

        public String name;
        public String expression;
        public String description;
        public String type;
        public String validExpression;
        public String noDataValue;
        public Integer spectralBandIndex;
        public Float spectralWavelength;
        public Float spectralBandwidth;
    }

    public static class Variable {

        public String name;
        public String type;
        public String value;
    }

    @TargetProduct
    private Product targetProduct;

    @SourceProducts
    private Product[] sourceProducts;

    @Parameter(alias = "targetBands", itemAlias = "targetBand",
               description = "List of descriptors defining the target bands.")
    private BandDescriptor[] targetBandDescriptors;
    @Parameter(alias = "variables", itemAlias = "variable",
               description = "List of variables which can be used within the expressions.")
    private Variable[] variables;

    private Map<Band, BandDescriptor> descriptorMap;

    public static BandMathsOp createBooleanExpressionBand(String expression, Product sourceProduct) {
        BandDescriptor[] bandDescriptors = new BandDescriptor[1];
        bandDescriptors[0] = new BandDescriptor();
        bandDescriptors[0].name = "band1";
        bandDescriptors[0].expression = expression;
        bandDescriptors[0].type = ProductData.TYPESTRING_INT8;

        BandMathsOp bandMathOp = new BandMathsOp();
        bandMathOp.targetBandDescriptors = bandDescriptors;
        bandMathOp.sourceProducts = new Product[]{sourceProduct};
        return bandMathOp;
    }

    @Override
    public void initialize() throws OperatorException {
        if (targetBandDescriptors == null || targetBandDescriptors.length == 0) {
            throw new OperatorException("No target bands specified.");
        }

        int width = sourceProducts[0].getSceneRasterWidth();
        int height = sourceProducts[0].getSceneRasterHeight();
        for (Product product : sourceProducts) {
            if (product.getSceneRasterWidth() != width ||
                product.getSceneRasterHeight() != height) {
                throw new OperatorException("Products must have the same raster dimension.");
            }
        }
        targetProduct = new Product(sourceProducts[0].getName() + "BandMath", "BandMath", width, height);

        descriptorMap = new HashMap<Band, BandDescriptor>(targetBandDescriptors.length);
        Namespace namespace = createNamespace();
        Parser verificationParser = new ParserImpl(namespace, true);
        for (BandDescriptor bandDescriptor : targetBandDescriptors) {
            createBand(bandDescriptor, verificationParser);
        }

        ProductUtils.copyMetadata(sourceProducts[0], targetProduct);
        ProductUtils.copyGeoCoding(sourceProducts[0], targetProduct);
    }

    @Override
    public void computeTile(Band band, Tile targetTile, ProgressMonitor pm) throws OperatorException {
        Rectangle rect = targetTile.getRectangle();
        Term term = createTerm(descriptorMap.get(band).expression);
        RasterDataSymbol[] refRasterDataSymbols = BandArithmetic.getRefRasterDataSymbols(term);

        for (RasterDataSymbol symbol : refRasterDataSymbols) {
            Tile tile = getSourceTile(symbol.getRaster(), rect, pm);
            if (tile.getRasterDataNode().isScalingApplied()) {
                ProductData dataBuffer = ProductData.createInstance(ProductData.TYPE_FLOAT32,
                                                                    tile.getWidth() * tile.getHeight());
                int dataBufferIndex = 0;
                for (int y = rect.y; y < rect.y + rect.height; y++) {
                    for (int x = rect.x; x < rect.x + rect.width; x++) {
                        dataBuffer.setElemFloatAt(dataBufferIndex, tile.getSampleFloat(x, y));
                        dataBufferIndex++;
                    }
                }
                symbol.setData(dataBuffer);
            } else {
                ProductData dataBuffer = tile.getRawSamples();
                symbol.setData(dataBuffer);
            }
        }

        final RasterDataEvalEnv env = new RasterDataEvalEnv(rect.x, rect.y, rect.width, rect.height);
        pm.beginTask("Evaluating expression", rect.height);
        try {
            int pixelIndex = 0;
            for (int y = rect.y; y < rect.y + rect.height; y++) {
                if (pm.isCanceled()) {
                    break;
                }
                for (int x = rect.x; x < rect.x + rect.width; x++) {
                    env.setElemIndex(pixelIndex);
                    targetTile.setSample(x, y, term.evalD(env));
                    pixelIndex++;
                }
                pm.worked(1);
            }
        } finally {
            pm.done();
        }
    }

    private void createBand(BandDescriptor bandDescriptor, Parser verificationParser) {
        if(StringUtils.isNullOrEmpty(bandDescriptor.name)) {
             throw new OperatorException("Missing band name.");
        }
        if(StringUtils.isNullOrEmpty(bandDescriptor.type)) {
             throw new OperatorException(String.format("Missing data type for band %s.", bandDescriptor.name));
        }

        Band band = targetProduct.addBand(bandDescriptor.name, ProductData.getType(bandDescriptor.type.toLowerCase()));
        if (StringUtils.isNotNullAndNotEmpty(bandDescriptor.description)) {
            band.setDescription(bandDescriptor.description);
        }
        if (StringUtils.isNotNullAndNotEmpty(bandDescriptor.validExpression)) {
            band.setValidPixelExpression(bandDescriptor.validExpression);
        }
        if (StringUtils.isNotNullAndNotEmpty(bandDescriptor.noDataValue)) {
            try {
                double parseDouble = Double.parseDouble(bandDescriptor.noDataValue);
                band.setNoDataValue(parseDouble);
                band.setNoDataValueUsed(true);
            } catch (NumberFormatException e) {
                throw new OperatorException("Bad value for NoDataValue given: " + bandDescriptor.noDataValue, e);
            }
        }
        if (bandDescriptor.spectralBandIndex != null) {
            band.setSpectralBandIndex(bandDescriptor.spectralBandIndex);
        }
        if (bandDescriptor.spectralWavelength != null) {
            band.setSpectralWavelength(bandDescriptor.spectralWavelength);
        }
        if (bandDescriptor.spectralBandwidth != null) {
            band.setSpectralBandwidth(bandDescriptor.spectralBandwidth);
        }
        descriptorMap.put(band, bandDescriptor);
        try {
            Term testTerm = verificationParser.parse(bandDescriptor.expression);
        } catch (ParseException e) {
            throw new OperatorException("Could not parse expression: " + bandDescriptor.expression, e);
        }
    }

    private Namespace createNamespace() {
        WritableNamespace namespace = BandArithmetic.createDefaultNamespace(sourceProducts, 0,
                                                                            new SourceProductPrefixProvider());
        if (variables != null) {
            for (Variable variable : variables) {
                if (ProductData.isFloatingPointType(ProductData.getType(variable.type))) {
                    Symbol symbol = SymbolFactory.createConstant(variable.name, Double.parseDouble(variable.value));
                    namespace.registerSymbol(symbol);
                } else if (ProductData.isIntType(ProductData.getType(variable.type))) {
                    Symbol symbol = SymbolFactory.createConstant(variable.name, Long.parseLong(variable.value));
                    namespace.registerSymbol(symbol);
                } else if ("boolean".equals(variable.type)) {
                    Symbol symbol = SymbolFactory.createConstant(variable.name, Boolean.parseBoolean(variable.value));
                    namespace.registerSymbol(symbol);
                }
            }
        }
        return namespace;
    }

    private Term createTerm(String expression) {
        Namespace namespace = createNamespace();
        final Term term;
        try {
            Parser parser = new ParserImpl(namespace, false);
            term = parser.parse(expression);
        } catch (ParseException e) {
            throw new OperatorException("Could not parse expression: " + expression, e);
        }
        return term;
    }

    public static class Spi extends OperatorSpi {

        public Spi() {
            super(BandMathsOp.class);
            setOperatorUI(BandArithmeticOpUI.class);
        }
    }

    private class SourceProductPrefixProvider implements ProductPrefixProvider {

        @Override
        public String getPrefix(Product product) {
            return "$" + getSourceProductId(product) + ".";
        }
    }
}
