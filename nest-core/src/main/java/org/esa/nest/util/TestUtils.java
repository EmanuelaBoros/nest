/*
 * Copyright (C) 2010 Array Systems Computing Inc. http://www.array.ca
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
package org.esa.nest.util;

import com.bc.ceres.core.ProgressMonitor;
import org.esa.beam.framework.dataio.*;
import org.esa.beam.framework.datamodel.*;
import org.esa.beam.framework.dataop.maptransf.Datum;
import org.esa.beam.framework.gpf.Operator;
import org.esa.beam.framework.gpf.OperatorSpi;
import org.esa.beam.gpf.operators.standard.WriteOp;
import org.esa.beam.util.PropertyMap;
import org.esa.beam.util.ProductUtils;
import org.esa.beam.dataio.dimap.DimapProductConstants;
import org.esa.nest.dataio.ReaderUtils;
import org.esa.nest.datamodel.AbstractMetadata;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;

/**
 * Utilities for Operator unit tests
 * In order to test the datasets at Array set the foloowing to true in the nest.config
 * nest.test.ReadersOnAllProducts=true nest.test.ProcessingOnAllProducts=true
 */
public class TestUtils {

    private static final PropertyMap testPreferences = createTestPreferences();

    private final static String contextID = ResourceUtils.getContextID();
    public final static String rootPathExpectedProducts = testPreferences.getPropertyString(contextID+".test.rootPathExpectedProducts");
    public final static String rootPathTerraSarX = testPreferences.getPropertyString(contextID+".test.rootPathTerraSarX");
    public final static String rootPathASAR= testPreferences.getPropertyString(contextID+".test.rootPathASAR");
    public final static String rootPathRadarsat2 = testPreferences.getPropertyString(contextID+".test.rootPathRadarsat2");
    public final static String rootPathRadarsat1 = testPreferences.getPropertyString(contextID+".test.rootPathRadarsat1");
    public final static String rootPathERS = testPreferences.getPropertyString(contextID+".test.rootPathERS");
    public final static String rootPathJERS = testPreferences.getPropertyString(contextID+".test.rootPathJERS");
    public final static String rootPathALOS = testPreferences.getPropertyString(contextID+".test.rootPathALOS");
    public final static String rootPathCosmoSkymed = testPreferences.getPropertyString(contextID+".test.rootPathCosmoSkymed");
    public final static String rootPathMixProducts = testPreferences.getPropertyString(contextID+".test.rootPathMixProducts");

    public final static int subsetX = Integer.parseInt(testPreferences.getPropertyString(contextID+".test.subsetX"));
    public final static int subsetY = Integer.parseInt(testPreferences.getPropertyString(contextID+".test.subsetY"));
    public final static int subsetWidth = Integer.parseInt(testPreferences.getPropertyString(contextID+".test.subsetWidth"));
    public final static int subsetHeight = Integer.parseInt(testPreferences.getPropertyString(contextID+".test.subsetHeight"));

    private static String[] nonValidExtensions = { "xsd", "xsl", "xls", "pdf", "txt", "doc", "ps", "db", "ief", "ord",
                                                   "tfw", "gif", "jpg", "jgw", "hdr", "self", "report", "raw", "tgz",
                                                   "log", "html", "htm", "png", "bmp", "ps", "aux", "ovr", "brs", "kml", "kmz",
                                                   "sav", "7z", "zip", "rrd", "lbl", "z", "gz", "exe", "bat", "sh", "rtf",
                                                   "prj", "dbf", "shx", "ace", "ace2", "tar"};
    private static String[] nonValidprefixes = { "led", "trl", "tra_", "nul", "lea", "dat", "img", "dfas", "dfdn", "lut",
                                                 "readme", "l1b_iif", "dor_vor", "imagery_", "browse" };

    private static final int maxIteration = Integer.parseInt(testPreferences.getPropertyString(contextID+".test.maxProductsPerRootFolder"));

    private static PropertyMap createTestPreferences() {
        final PropertyMap prefs = new PropertyMap();
        try {
            prefs.load(ResourceUtils.findConfigFile(ResourceUtils.getContextID()+".config"));
        } catch(IOException e) {
            System.out.println("Unable to load test preferences "+e.getMessage());
        }
        return prefs;
    }

    public static int getMaxIterations() {
        return maxIteration;
    }

    public static boolean canTestReadersOnAllProducts() {
        final String testAllProducts = testPreferences.getPropertyString(contextID+".test.ReadersOnAllProducts");
        return testAllProducts != null && testAllProducts.equalsIgnoreCase("true");
    }

    public static boolean canTestProcessingOnAllProducts() {
        final String testAllProducts = testPreferences.getPropertyString(contextID+".test.ProcessingOnAllProducts");
        return testAllProducts != null && testAllProducts.equalsIgnoreCase("true");
    }

    public static Product createProduct(final String type, final int w, final int h) {
        Product product = new Product("name", type, w, h);

        product.setStartTime(AbstractMetadata.parseUTC("10-MAY-2008 20:30:46.890683"));
        product.setEndTime(AbstractMetadata.parseUTC("10-MAY-2008 20:35:46.890683"));
        product.setDescription("description");

        addGeoCoding(product);

        AbstractMetadata.addAbstractedMetadataHeader(product.getMetadataRoot());

        return product;
    }

    private static void addGeoCoding(final Product product) {

        TiePointGrid latGrid = new TiePointGrid("lat", 2, 2, 0.5f, 0.5f,
                product.getSceneRasterWidth(), product.getSceneRasterHeight(),
                      new float[]{10.0f, 10.0f, 5.0f, 5.0f});
        TiePointGrid lonGrid = new TiePointGrid("lon", 2, 2, 0.5f, 0.5f,
                product.getSceneRasterWidth(), product.getSceneRasterHeight(),
                      new float[]{10.0f, 10.0f, 5.0f, 5.0f},
                      TiePointGrid.DISCONT_AT_360);
        TiePointGeoCoding tpGeoCoding = new TiePointGeoCoding(latGrid, lonGrid, Datum.WGS_84);

        product.addTiePointGrid(latGrid);
        product.addTiePointGrid(lonGrid);
        product.setGeoCoding(tpGeoCoding);
    }

    public static void verifyProduct(final Product product, final boolean verifyTimes,
                                     final boolean verifyGeoCoding) throws Exception {
        ReaderUtils.verifyProduct(product, verifyTimes, verifyGeoCoding);
    }

    public static void attributeEquals(final MetadataElement elem, final String name,
                                       final double trueValue) throws Exception {
        double val = elem.getAttributeDouble(name, 0);
        if(Double.compare(val, trueValue) != 0) {
            if(Float.compare((float)val, (float)trueValue) != 0)
                throwErr(name + " is " + val + ", expecting " + trueValue);
        }
    }

    public static void attributeEquals(final MetadataElement elem, String name,
                                       final String trueValue) throws Exception {
        String val = elem.getAttributeString(name, "");
        if(!val.equals(trueValue))
            throwErr(name + " is " + val + ", expecting " + trueValue);
    }

    private static void compareMetadata(final Product testProduct, final Product expectedProduct,
                                        final String[] excemptionList) throws Exception {
        final MetadataElement testAbsRoot = AbstractMetadata.getAbstractedMetadata(testProduct);
        if(testAbsRoot == null)
            throwErr("Metadata is null");
        final MetadataElement expectedAbsRoot = AbstractMetadata.getAbstractedMetadata(expectedProduct);
        if(expectedAbsRoot == null)
            throwErr("Metadata is null");

        if(excemptionList != null) {
            Arrays.sort(excemptionList);
        }

        final MetadataAttribute[] attribList = expectedAbsRoot.getAttributes();
        for(MetadataAttribute expectedAttrib : attribList) {
            if(excemptionList != null && Arrays.binarySearch(excemptionList, expectedAttrib.getName()) >= 0)
                continue;

            final MetadataAttribute result = testAbsRoot.getAttribute(expectedAttrib.getName());
            if(result == null) {
                throwErr("Metadata attribute "+expectedAttrib.getName()+" is missing");
            }
            if(!result.getData().equalElems(expectedAttrib.getData())) {
                if(expectedAttrib.getData().toString().trim().equalsIgnoreCase(result.getData().toString().trim())) {

                } else {
                    throwErr("Metadata attribute "+expectedAttrib.getName()+" expecting "+expectedAttrib.getData().toString()
                        +" got "+ result.getData().toString());
                }
            }
        }
    }

    public static void compareProducts(final Operator op, final Product targetProduct,
                                       final String expectedPath, final String[] excemptionList) throws Exception {

        final Band targetBand = targetProduct.getBandAt(0);
        if(targetBand == null)
            throwErr("targetBand at 0 is null");

        // readPixels: execute computeTiles()
        final float[] floatValues = new float[10000];
        targetBand.readPixels(100, 101, 100, 99, floatValues, ProgressMonitor.NULL);

        // compare with expected outputs:
        final File expectedFile = new File(expectedPath);
        if(!expectedFile.exists()) {
            throwErr("Expected file not found "+expectedFile.toString());
        }

        final ProductReader reader2 = ProductIO.getProductReaderForFile(expectedFile);

        final Product expectedProduct = reader2.readProductNodes(expectedFile, null);
        final Band expectedBand = expectedProduct.getBandAt(0);

        final float[] expectedValues = new float[10000];
        expectedBand.readPixels(100, 101, 100, 99, expectedValues, ProgressMonitor.NULL);
        if(!Arrays.equals(floatValues, expectedValues))
                throwErr("Pixels are different");

        // compare updated metadata
        compareMetadata(targetProduct, expectedProduct, excemptionList);
    }

    public static void executeOperator(final Operator op) throws Exception {
        // get targetProduct: execute initialize()
        final Product targetProduct = op.getTargetProduct();
        TestUtils.verifyProduct(targetProduct, false, !isAlos(targetProduct));

        final Band targetBand = targetProduct.getBandAt(0);
        if(targetBand == null)
            throwErr("targetBand at 0 is null");

        final int bandWidth = targetBand.getSceneRasterWidth();
        final int bandHeight = targetBand.getSceneRasterHeight();

        // readPixels: execute computeTiles()
        final float[] floatValues = new float[10000];
        targetBand.readPixels(within(subsetX, bandWidth),
                              within(subsetY, bandHeight),
                              within(subsetWidth, bandWidth),
                              within(subsetHeight, bandHeight), 
                              floatValues, ProgressMonitor.NULL);
    }

    public static Product createSubsetProduct(final Product sourceProduct) throws IOException {
        final int bandWidth = sourceProduct.getSceneRasterWidth();
        final int bandHeight = sourceProduct.getSceneRasterHeight();

        final ProductSubsetBuilder subsetReader = new ProductSubsetBuilder();
        final ProductSubsetDef subsetDef = new ProductSubsetDef();

        subsetDef.addNodeNames(sourceProduct.getTiePointGridNames());

        final String bandName = ProductUtils.findSuitableQuicklookBandName(sourceProduct);
        subsetDef.addNodeNames(new String[] { bandName } );
        subsetDef.setRegion(within(subsetX, bandWidth),
                            within(subsetY, bandHeight),
                            within(subsetWidth, bandWidth),
                            within(subsetHeight, bandHeight));
        subsetDef.setIgnoreMetadata(false);
        subsetDef.setTreatVirtualBandsAsRealBands(true);

        final Product subsetProduct = subsetReader.readProductNodes(sourceProduct, subsetDef);
        final File tmpFile = new File(ResourceUtils.getApplicationUserTempDataDir(), "tmp_subset.dim");
        WriteOp.writeProduct(subsetProduct, tmpFile, DimapProductConstants.DIMAP_FORMAT_NAME, ProgressMonitor.NULL);           

        return ProductIO.readProduct(tmpFile);
    }

    private static int within(final int val, final int max) {
        return Math.max(0, Math.min(val, max));
    }

    private static boolean isAlos(Product prod) {
        final MetadataElement absRoot = AbstractMetadata.getAbstractedMetadata(prod);
        if(absRoot != null) {
            return absRoot.getAttributeString(AbstractMetadata.MISSION).contains("ALOS");
        }
        return false;
    }

    public static int recurseProcessFolder(final OperatorSpi spi, final File folder, int iterations,
                                            final String[] productTypeExemptions,
                                            final String[] exceptionExemptions) throws Exception {
        for(File file : folder.listFiles()) {
            if(maxIteration > 0 && iterations >= maxIteration)
                break;

            if(file.isDirectory()) {
                if(!file.getName().contains("skipTest")) {
                    iterations = recurseProcessFolder(spi, file, iterations, productTypeExemptions, exceptionExemptions);
                }
            } else {
                try {
                    if(isNotProduct(file))
                        continue;
                    final ProductReader reader = ProductIO.getProductReaderForFile(file);
                    if(reader != null) {
                        final Product sourceProduct = reader.readProductNodes(file, null);
                        if(contains(sourceProduct.getProductType(), productTypeExemptions))
                            continue;

                        final Operator op = spi.createOperator();
                        op.setSourceProduct(sourceProduct);

                        System.out.println(spi.getOperatorAlias()+" Processing "+ file.toString());
                        TestUtils.executeOperator(op);

                        ++iterations;
                    } else {
                        System.out.println(file.getName() + " is non valid");
                    }
                } catch(Exception e) {
                    boolean ok = false;
                    if(exceptionExemptions != null) {
                        for(String excemption : exceptionExemptions) {
                            if(e.getMessage().contains(excemption)) {
                                ok = true;
                                System.out.println("Excemption for "+e.getMessage());
                                break;
                            }
                        }
                    }
                    if(!ok) {
                        //System.out.println("Failed to process "+ file.toString());
                        throw e;
                    }
                }
            }
        }
        return iterations;
    }

    public static boolean isNotProduct(final File file) {
        final String name = file.getName().toLowerCase();
        for(String ext : nonValidExtensions) {
            if(name.endsWith(ext))
                return true;
        }
        for(String pre : nonValidprefixes) {
            if(name.startsWith(pre))
                return true;
        }
        return false;
    }

    public static boolean contains(final String value, final String[] exemptions) {
        if(exemptions != null) {
            for(String type : exemptions) {
                if(value.contains(type))
                    return true;
            }
        }
        return false;
    }

    /**
     * Processes all products in a folder
     * @param spi the OperatorSpi to create the operator
     * @param folderPath the path to recurse through
     * @param productTypeExemptions product types to ignore
     * @param exceptionExemptions exceptions that are ok and can be ignored for the test
     * @throws Exception general exception
     */
    public static void testProcessAllInPath(final OperatorSpi spi, final String folderPath,
                                            final String[] productTypeExemptions,
                                            final String[] exceptionExemptions) throws Exception
    {
        final File folder = new File(folderPath);
        if(!folder.exists()) return;

        if(canTestProcessingOnAllProducts()) {
            int iterations = 0;
            recurseProcessFolder(spi, folder, iterations, productTypeExemptions, exceptionExemptions);
        }
    }

    /**
     * Processes all products in a folder
     * @param processor the RecursiveProcessor to create the graph
     * @param folderPath the path to recurse through
     * @param productTypeExemptions product types to ignore
     * @param exceptionExemptions exceptions that are ok and can be ignored for the test
     * @throws Exception general exception
     */
    public static void testProcessAllInPath(final RecursiveProcessor processor, final String folderPath,
                                            final String[] productTypeExemptions,
                                            final String[] exceptionExemptions) throws Exception
    {
        final File folder = new File(folderPath);
        if(!folder.exists()) return;

        if(canTestProcessingOnAllProducts()) {
            int iterations = 0;
            processor.recurseProcessFolder(folder, iterations, productTypeExemptions, exceptionExemptions);
        }
    }

    public static void recurseReadFolder(final File folder,
                                         final ProductReaderPlugIn readerPlugin,
                                         final ProductReader reader,
                                         final String[] productTypeExemptions,
                                         final String[] exceptionExemptions) throws Exception {
        for(File file : folder.listFiles()) {
            if(file.isDirectory()) {
                if(!file.getName().contains("skipTest")) {
                    recurseReadFolder(file, readerPlugin, reader, productTypeExemptions, exceptionExemptions);
                }
            } else if(readerPlugin.getDecodeQualification(file) == DecodeQualification.INTENDED) {

                try {
                    System.out.println("Reading "+ file.toString());

                    final Product product = reader.readProductNodes(file, null);
                    if(contains(product.getProductType(), productTypeExemptions))
                            continue;
                    ReaderUtils.verifyProduct(product, true);
                } catch(Exception e) {
                    boolean ok = false;
                    if(exceptionExemptions != null) {
                        for(String excemption : exceptionExemptions) {
                            if(e.getMessage().contains(excemption)) {
                                ok = true;
                                System.out.println("Excemption for "+e.getMessage());
                                break;
                            }
                        }
                    }
                    if(!ok) {
                        System.out.println("Failed to read "+ file.toString());
                        throw e;
                    }
                }
            }
        }
    }

    private static void throwErr(final String description) throws Exception {
        throw new Exception(description);
    }
}