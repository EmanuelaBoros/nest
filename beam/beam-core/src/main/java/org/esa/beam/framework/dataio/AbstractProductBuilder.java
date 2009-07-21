/*
 * $id$
 *
 * Copyright (c) 2003 Brockmann Consult GmbH. All right reserved.
 * http://www.brockmann-consult.de
 */
package org.esa.beam.framework.dataio;

import org.esa.beam.framework.datamodel.*;
import org.esa.beam.util.Guardian;
import org.esa.beam.util.ProductUtils;

import java.io.IOException;
import java.util.Hashtable;
import java.util.Map;

public abstract class AbstractProductBuilder extends AbstractProductReader {

    protected Product sourceProduct;
    protected boolean sourceProductOwner;
    protected int sceneRasterWidth;
    protected int sceneRasterHeight;
    protected String newProductName;
    protected String newProductDesc;
    protected Map<Band, RasterDataNode> bandMap;

    public AbstractProductBuilder(final boolean sourceProductOwner) {
        super(null);
        bandMap = new Hashtable<Band, RasterDataNode>(16);
        this.sourceProductOwner = sourceProductOwner;
    }

    public Product getSourceProduct() {
        return sourceProduct;
    }

    public boolean isSourceProductOwner() {
        return sourceProductOwner;
    }

    public void setNewProductDesc(String newProductDesc) {
        this.newProductDesc = newProductDesc;
    }

    public void setNewProductName(String newProductName) {
        this.newProductName = newProductName;
    }

    public int getSceneRasterWidth() {
        return sceneRasterWidth;
    }

    public int getSceneRasterHeight() {
        return sceneRasterHeight;
    }

    protected Product readProductNodes(Product sourceProduct, ProductSubsetDef subsetDef, String name,
                                       String desc) throws IOException {
        Guardian.assertNotNull("sourceProduct", sourceProduct);
        setNewProductName(name != null ? name : sourceProduct.getName());
        setNewProductDesc(desc != null ? desc : sourceProduct.getDescription());
        final Product product = readProductNodes(sourceProduct, subsetDef);
        if (sourceProduct.getQuicklookBandName() != null
                && product.getQuicklookBandName() == null
                && product.containsBand(sourceProduct.getQuicklookBandName())) {
            product.setQuicklookBandName(sourceProduct.getQuicklookBandName());
        }
        product.setModified(true);
        return product;
    }

    @Override
    protected abstract Product readProductNodesImpl() throws IOException;

    /**
     * Closes the access to all currently opened resources such as file input streams and all resources of this children
     * directly owned by this reader. Its primary use is to allow the garbage collector to perform a vanilla job.
     * <p/>
     * <p>This method should be called only if it is for sure that this object instance will never be used again. The
     * results of referencing an instance of this class after a call to <code>close()</code> are undefined.
     * <p/>
     * <p>Overrides of this method should always call <code>super.close();</code> after disposing this instance.
     *
     * @throws IOException if an I/O error occurs
     */
    @Override
    public void close() throws IOException {
        disposeBandMap();
        if (sourceProductOwner && sourceProduct != null) {
            sourceProduct.dispose();
        }
        sourceProduct = null;
        super.close();
    }

    protected void addFlagCodingsToProduct(Product product) {
        final ProductNodeGroup<FlagCoding> flagCodingGroup = getSourceProduct().getFlagCodingGroup();
        for (int i = 0; i < flagCodingGroup.getNodeCount(); i++) {
            FlagCoding sourceFlagCoding = flagCodingGroup.get(i);
            FlagCoding destFlagCoding = new FlagCoding(sourceFlagCoding.getName());
            destFlagCoding.setDescription(sourceFlagCoding.getDescription());
            cloneFlags(sourceFlagCoding, destFlagCoding);
            product.getFlagCodingGroup().add(destFlagCoding);
        }
    }

    protected void addIndexCodingsToProduct(Product product) {
        final ProductNodeGroup<IndexCoding> indexCodingGroup = getSourceProduct().getIndexCodingGroup();
        for (int i = 0; i < indexCodingGroup.getNodeCount(); i++) {
            IndexCoding sourceIndexCoding = indexCodingGroup.get(i);
            IndexCoding destIndexCoding = new IndexCoding(sourceIndexCoding.getName());
            destIndexCoding.setDescription(sourceIndexCoding.getDescription());
            cloneIndexes(sourceIndexCoding, destIndexCoding);
            product.getIndexCodingGroup().add(destIndexCoding);
        }
    }

    protected static void addAttribString(String name, String value, MetadataElement subsetElem) {
        final ProductData data = ProductData.createInstance(value);
        subsetElem.addAttributeFast(new MetadataAttribute(name, data, true));
    }

    protected void addBitmaskDefsToProduct(Product product) {
        ProductUtils.copyBitmaskDefsAndOverlays(getSourceProduct(), product);
    }

    protected void addBitmaskOverlayInfosToBandAndTiePointGrids(final Product product) {
        copyBitmaskOverlayInfo(getSourceProduct().getBands(), product);
        copyBitmaskOverlayInfo(getSourceProduct().getTiePointGrids(), product);
    }

    private static void copyBitmaskOverlayInfo(final RasterDataNode[] sourceNodes, final Product product) {
        for (final RasterDataNode sourceNode : sourceNodes) {
            final RasterDataNode destNode = product.getRasterDataNode(sourceNode.getName());
            if (destNode != null) {
                final BitmaskOverlayInfo bitmaskOverlayInfo = sourceNode.getBitmaskOverlayInfo();
                if (bitmaskOverlayInfo != null) {
                    final BitmaskOverlayInfo info = new BitmaskOverlayInfo();
                    final BitmaskDef[] bitmaskDefs = bitmaskOverlayInfo.getBitmaskDefs();
                    for (BitmaskDef bitmaskDef : bitmaskDefs) {
                        info.addBitmaskDef(product.getBitmaskDef(bitmaskDef.getName()));
                    }
                    destNode.setBitmaskOverlayInfo(info);
                }
            }
        }
    }

    protected void cloneFlags(FlagCoding sourceFlagCoding, FlagCoding destFlagCoding) {
        cloneMetadataElementsAndAttributes(sourceFlagCoding, destFlagCoding, 1);
    }

    protected void cloneIndexes(IndexCoding sourceFlagCoding, IndexCoding destFlagCoding) {
        cloneMetadataElementsAndAttributes(sourceFlagCoding, destFlagCoding, 1);
    }

    protected void addMetadataToProduct(Product product) {
        cloneMetadataElementsAndAttributes(getSourceProduct().getMetadataRoot(), product.getMetadataRoot(), 1);
    }

    protected void cloneMetadataElementsAndAttributes(MetadataElement sourceRoot, MetadataElement destRoot, int level) {
        cloneMetadataElements(sourceRoot, destRoot, level);
        cloneMetadataAttributes(sourceRoot, destRoot);
    }

    protected void cloneMetadataElements(MetadataElement sourceRoot, MetadataElement destRoot, int level) {
        for (int i = 0; i < sourceRoot.getNumElements(); i++) {
            MetadataElement sourceElement = sourceRoot.getElementAt(i);
            if (level > 0 || isNodeAccepted(sourceElement.getName())) {
                MetadataElement element = new MetadataElement(sourceElement.getName());
                element.setDescription(sourceElement.getDescription());
                destRoot.addElement(element);
                cloneMetadataElementsAndAttributes(sourceElement, element, level + 1);
            }
        }
    }

    protected void cloneMetadataAttributes(MetadataElement sourceRoot, MetadataElement destRoot) {
        for (int i = 0; i < sourceRoot.getNumAttributes(); i++) {
            final MetadataAttribute sourceAttribute = sourceRoot.getAttributeAt(i);
            destRoot.addAttributeFast(sourceAttribute.createDeepClone());
        }
    }

    @Override
    protected boolean isInstanceOfValidInputType(Object input) {
        return input instanceof Product;
    }

    protected void disposeBandMap() {
        bandMap.clear();
    }


}
