/*
 * $Id: MergeBandsOp.java,v 1.1 2008-09-25 21:26:07 lveci Exp $
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
package org.esa.nest.gpf;

import com.bc.ceres.core.ProgressMonitor;
import org.esa.beam.framework.datamodel.Band;
import org.esa.beam.framework.datamodel.FlagCoding;
import org.esa.beam.framework.datamodel.IndexCoding;
import org.esa.beam.framework.datamodel.Product;
import org.esa.beam.framework.gpf.Operator;
import org.esa.beam.framework.gpf.OperatorException;
import org.esa.beam.framework.gpf.OperatorSpi;
import org.esa.beam.framework.gpf.Tile;
import org.esa.beam.framework.gpf.operators.common.MergeOp;
import org.esa.beam.framework.gpf.annotations.OperatorMetadata;
import org.esa.beam.framework.gpf.annotations.Parameter;
import org.esa.beam.framework.gpf.annotations.TargetProduct;
import org.esa.beam.framework.gpf.annotations.SourceProducts;
import org.esa.beam.util.ProductUtils;
import org.esa.beam.util.StringUtils;

import java.awt.image.RenderedImage;
import java.util.HashSet;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

@OperatorMetadata(alias = "Merge-Bands",
        description = "Merges an arbitrary number of source bands into the target product.")
public class MergeBandsOp extends Operator {

    @SourceProducts
    private Product[] sourceProducts;


    @Parameter(defaultValue = "mergedProduct", description = "The name of the target product.")
    private String productName;

    @TargetProduct
    private Product targetProduct;


    @Override
    public void initialize() throws OperatorException {

        Product srcProduct = sourceProducts[0];
        final int sceneRasterWidth = srcProduct.getSceneRasterWidth();
        final int sceneRasterHeight = srcProduct.getSceneRasterHeight();
        targetProduct = new Product(productName, srcProduct.getProductType(),
                sceneRasterWidth, sceneRasterHeight);

        copyGeoCoding(srcProduct, targetProduct);

     /*
        Set<Product> allSrcProducts = new HashSet<Product>();
        for (BandDesc bandDesc : bands) {
            Product srcProduct = getSourceProduct(bandDesc.product);
            if (StringUtils.isNotNullAndNotEmpty(bandDesc.name)) {
                if (StringUtils.isNotNullAndNotEmpty(bandDesc.newName)) {
                    copyBandWithFeatures(srcProduct, targetProduct, bandDesc.name, bandDesc.newName);
                } else {
                    copyBandWithFeatures(srcProduct, targetProduct, bandDesc.name);
                }
                allSrcProducts.add(srcProduct);
            } else if (StringUtils.isNotNullAndNotEmpty(bandDesc.nameExp)) {
                Pattern pattern = Pattern.compile(bandDesc.nameExp);
                for (String bandName : srcProduct.getBandNames()) {
                    Matcher matcher = pattern.matcher(bandName);
                    if (matcher.matches()) {
                        copyBandWithFeatures(srcProduct, targetProduct, bandName);
                        allSrcProducts.add(srcProduct);
                    }
                }
            }
        }

        for (Product srcProduct : allSrcProducts) {
            ProductUtils.copyBitmaskDefsAndOverlays(srcProduct, targetProduct);
        }   */
    }

    /*
     * Copies the tie point data, geocoding and the start and stop time.
     */
    private static void copyGeoCoding(Product sourceProduct,
                                      Product destinationProduct) {
        // copy all tie point grids to output product
        ProductUtils.copyTiePointGrids(sourceProduct, destinationProduct);
        // copy geo-coding to the output product
        ProductUtils.copyGeoCoding(sourceProduct, destinationProduct);
        destinationProduct.setStartTime(sourceProduct.getStartTime());
        destinationProduct.setEndTime(sourceProduct.getEndTime());
    }

    private void copyBandWithFeatures(Product srcProduct, Product outputProduct, String oldBandName, String newBandName) {
        Band destBand = copyBandWithFeatures(srcProduct, outputProduct, oldBandName);
        destBand.setName(newBandName);
    }

    private Band copyBandWithFeatures(Product srcProduct, Product outputProduct, String bandName) {
        Band destBand = ProductUtils.copyBand(bandName, srcProduct, outputProduct);
        Band srcBand = srcProduct.getBand(bandName);
        destBand.setSourceImage(srcBand.getSourceImage());
        if (srcBand.getFlagCoding() != null) {
            FlagCoding srcFlagCoding = srcBand.getFlagCoding();
            if (!outputProduct.getFlagCodingGroup().contains(srcFlagCoding.getName())) {
                ProductUtils.copyFlagCoding(srcFlagCoding, outputProduct);
            }
            destBand.setSampleCoding(outputProduct.getFlagCodingGroup().get(srcFlagCoding.getName()));
        }
        if (srcBand.getIndexCoding() != null) {
            IndexCoding srcIndexCoding = srcBand.getIndexCoding();
            if (!outputProduct.getIndexCodingGroup().contains(srcIndexCoding.getName())) {
                ProductUtils.copyIndexCoding(srcIndexCoding, outputProduct);
            }
            destBand.setSampleCoding(outputProduct.getIndexCodingGroup().get(srcIndexCoding.getName()));
        }
        return destBand;
    }

    @Override
    public void computeTile(Band band, Tile targetTile, ProgressMonitor pm) throws OperatorException {
        getLogger().warning("Wrongly configured ProductMerger operator. Tiles should not be requested.");
    }


    public static class Spi extends OperatorSpi {
        public Spi() {
            super(MergeBandsOp.class);
            setOperatorUI(MergeBandsOpUI.class);
        }
    }
}