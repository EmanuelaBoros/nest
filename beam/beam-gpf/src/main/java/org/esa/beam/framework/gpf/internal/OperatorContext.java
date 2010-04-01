/*
 * $Id: OperatorContext.java,v 1.14 2010-04-01 19:42:54 lveci Exp $
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
package org.esa.beam.framework.gpf.internal;

import com.bc.ceres.binding.ConversionException;
import com.bc.ceres.binding.Property;
import com.bc.ceres.binding.PropertyContainer;
import com.bc.ceres.binding.PropertyDescriptor;
import com.bc.ceres.binding.PropertyDescriptorFactory;
import com.bc.ceres.binding.ValidationException;
import com.bc.ceres.binding.ValueSet;
import com.bc.ceres.binding.dom.DefaultDomConverter;
import com.bc.ceres.binding.dom.DomElement;
import com.bc.ceres.binding.dom.Xpp3DomElement;
import com.bc.ceres.core.Assert;
import com.bc.ceres.core.ProgressMonitor;
import com.bc.ceres.glevel.MultiLevelImage;
import org.esa.beam.framework.dataio.ProductReader;
import org.esa.beam.framework.datamodel.Band;
import org.esa.beam.framework.datamodel.MetadataAttribute;
import org.esa.beam.framework.datamodel.MetadataElement;
import org.esa.beam.framework.datamodel.Product;
import org.esa.beam.framework.datamodel.ProductData;
import org.esa.beam.framework.datamodel.RasterDataNode;
import org.esa.beam.framework.gpf.GPF;
import org.esa.beam.framework.gpf.Operator;
import org.esa.beam.framework.gpf.OperatorException;
import org.esa.beam.framework.gpf.OperatorSpi;
import org.esa.beam.framework.gpf.Tile;
import org.esa.beam.framework.gpf.annotations.ParameterDescriptorFactory;
import org.esa.beam.framework.gpf.annotations.SourceProduct;
import org.esa.beam.framework.gpf.annotations.SourceProducts;
import org.esa.beam.framework.gpf.annotations.TargetProduct;
import org.esa.beam.framework.gpf.annotations.TargetProperty;
import org.esa.beam.framework.gpf.graph.GraphOp;
import org.esa.beam.framework.gpf.internal.OperatorConfiguration.Reference;
import org.esa.beam.util.jai.JAIUtils;
import org.esa.beam.util.logging.BeamLogManager;

import javax.media.jai.JAI;
import java.awt.Dimension;
import java.awt.Rectangle;
import java.awt.RenderingHints;
import java.awt.image.Raster;
import java.awt.image.RenderedImage;
import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.lang.reflect.Modifier;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.logging.Logger;
import java.util.regex.Pattern;

/**
 * The context in which operators are executed.
 *
 * @author Norman Fomferra
 * @since 4.1
 */
public class OperatorContext {

    private static final String OPERATION_CANCELED_MESSAGE = "Operation canceled.";
    private String id;
    private Product targetProduct;
    private OperatorSpi operatorSpi;
    private Operator operator;
    private boolean computeTileMethodUsable;
    private boolean computeTileStackMethodUsable;
    private boolean passThrough;
    private List<Product> sourceProductList;
    private Map<String, Object> parameters;
    private Map<String, Product> sourceProductMap;
    private Map<String, Object> targetPropertyMap;
    private Map<Band, OperatorImage> targetImageMap;
    private OperatorConfiguration configuration;
    private Logger logger;
    private boolean disposed;
    private PropertyContainer propertyContainer;
    private RenderingHints renderingHints;
    private PerformanceMetric performanceMetric;
    private boolean initialising;
    private boolean requiresAllBands;


    public OperatorContext(Operator operator) {
        this.operator = operator;
        this.computeTileMethodUsable = canOperatorComputeTile(operator.getClass());
        this.computeTileStackMethodUsable = canOperatorComputeTileStack(operator.getClass());
        this.sourceProductList = new ArrayList<Product>(3);
        this.sourceProductMap = new HashMap<String, Product>(3);
        this.targetPropertyMap = new HashMap<String, Object>(3);
        this.logger = BeamLogManager.getSystemLogger();
        // Note: All OperatorImages in this context share the same TILE_CACHE_METRIC instance!
        this.performanceMetric = new PerformanceMetric();
        this.renderingHints = new RenderingHints(JAI.KEY_TILE_CACHE_METRIC, performanceMetric);
    }

    public String getId() {
        if (id == null) {
            id = getOperatorSpi().getOperatorAlias() + "$" + Long.toHexString(System.currentTimeMillis()).toUpperCase();
        }
        return id;
    }

    public void setId(String id) {
        Assert.notNull(id, "id");
        this.id = id;
    }

    public Logger getLogger() {
        return logger;
    }

    public void setLogger(Logger logger) {
        Assert.notNull(logger, "logger");
        this.logger = logger;
    }

    public PerformanceMetric getPerformanceMetric() {
        return performanceMetric;
    }

    public Product getSourceProduct(String id) {
        return sourceProductMap.get(id);
    }

    public void setSourceProduct(String id, Product product) {
        if (product != null) {
            if (!sourceProductList.contains(product)) {
                sourceProductList.add(product);
            }
            sourceProductMap.put(id, product);
        }
    }

    public Product[] getSourceProducts() {
        return sourceProductList.toArray(new Product[sourceProductList.size()]);
    }

    public void setSourceProducts(Product[] products) {
        sourceProductList.clear();
        sourceProductMap.clear();
        for (int i = 0; i < products.length; i++) {
            Product product = products[i];
            setSourceProduct(GPF.SOURCE_PRODUCT_FIELD_NAME + (i + 1), product);
        }
    }

    public void setSourceProducts(Map<String, Product> sourceProducts) {
        sourceProductList.clear();
        sourceProductMap.clear();
        Set<Map.Entry<String, Product>> entries = sourceProducts.entrySet();
        for (Map.Entry<String, Product> entry : entries) {
            setSourceProduct(entry.getKey(), entry.getValue());
        }
    }

    public String getSourceProductId(Product product) {
        Set<Map.Entry<String, Product>> entrySet = sourceProductMap.entrySet();
        for (Map.Entry<String, Product> entry : entrySet) {
            if (entry.getValue() == product) {
                return entry.getKey();
            }
        }
        return null;
    }

    public Product getTargetProduct() throws OperatorException {
        if (targetProduct == null) {
            initializeOperator();
        }
        return targetProduct;
    }

    public void setTargetProduct(Product targetProduct) {
        Assert.notNull(targetProduct, "targetProduct");
        this.targetProduct = targetProduct;
    }

    public Object getTargetProperty(String name) {
        getTargetProduct();

        return targetPropertyMap.get(name);
    }

    public boolean isPassThrough() {
        return passThrough;
    }

    public static void checkForCancelation(ProgressMonitor pm) throws OperatorException {
        if (pm.isCanceled()) {
            throw new OperatorException(OPERATION_CANCELED_MESSAGE);
        }
    }

    public OperatorSpi getOperatorSpi() {
        if (operatorSpi == null) {
            // create anonymous SPI
            operatorSpi = new OperatorSpi(operator.getClass()) {
            };
        }
        return operatorSpi;
    }

    public void setOperatorSpi(OperatorSpi operatorSpi) {
        this.operatorSpi = operatorSpi;
    }

    public Operator getOperator() {
        return operator;
    }


    public Object getParameter(String name) {
        Assert.notNull(name, "name");
        if (parameters == null) {
            return null;
        }
        return parameters.get(name);
    }

    public void setParameter(String name, Object value) {
        Assert.notNull(name, "name");
        if (value != null) {
            if (parameters == null) {
                parameters = new HashMap<String, Object>();
            }
            parameters.put(name, value);
        } else if (parameters != null) {
            parameters.remove(name);
        }
    }

    public Map<String, Object> getParameters() {
        return parameters;
    }

    public void setParameters(Map<String, Object> parameters) {
        this.parameters = new HashMap<String, Object>(parameters);
    }

    public RenderingHints getRenderingHints() {
        return renderingHints;
    }

    public void addRenderingHints(RenderingHints renderingHints) {
        this.renderingHints.add(renderingHints);
    }

    public void setConfiguration(OperatorConfiguration opConfiguration) {
        this.configuration = opConfiguration;
    }

    public boolean isInitialized() {
        return targetProduct != null;
    }

    public boolean isComputeTileMethodUsable() {
        return computeTileMethodUsable;
    }

    public boolean isComputeTileStackMethodUsable() {
        return computeTileStackMethodUsable;
    }

    public void setComputeTileMethodUsable(boolean computeTileMethodUsable) {
        this.computeTileMethodUsable = computeTileMethodUsable;
    }

    public void setComputeTileStackMethodUsable(boolean computeTileStackMethodUsable) {
        this.computeTileStackMethodUsable = computeTileStackMethodUsable;
    }

    public static Tile getSourceTile(RasterDataNode rasterDataNode, Rectangle rectangle, ProgressMonitor pm) {
        RenderedImage image = rasterDataNode.getSourceImage();
        ProgressMonitor oldPm = OperatorImage.setProgressMonitor(image, pm);
        try {
            /////////////////////////////////////////////////////////////////////
            //
            // Note: GPF pull-processing is triggered here!
            //
            Raster awtRaster = image.getData(rectangle); // Note: copyData is NOT faster!
            //
            /////////////////////////////////////////////////////////////////////
            return new TileImpl(rasterDataNode, awtRaster);
        } catch(Exception e) {
            System.out.println(e.toString());
            return null;    
        } finally {
            OperatorImage.setProgressMonitor(image, oldPm);
        }
    }

    public OperatorImage getTargetImage(Band band) {
        return targetImageMap.get(band);
    }

    public boolean isDisposed() {
        return disposed;
    }

    public void dispose() {
        if (!disposed) {
            disposed = true;
            parameters = null;
            configuration = null;
            sourceProductMap.clear();
            sourceProductList.clear();
            Collection<OperatorImage> operatorImages = targetImageMap.values();
            for (OperatorImage image : operatorImages) {
                image.dispose();
            }
            targetImageMap.clear();
            operator.dispose();
        }
    }

    private static boolean canOperatorComputeTile(Class<? extends Operator> aClass) {
        return implementsMethod(aClass, "computeTile",
                                new Class[]{
                                        Band.class,
                                        Tile.class,
                                        ProgressMonitor.class
                                });
    }

    private static boolean canOperatorComputeTileStack(Class<? extends Operator> aClass) {
        return implementsMethod(aClass, "computeTileStack",
                                new Class[]{
                                        Map.class,
                                        Rectangle.class,
                                        ProgressMonitor.class
                                });
    }

    private boolean operatorMustComputeTileStack() {
        return isComputeTileStackMethodUsable() && !isComputeTileMethodUsable();
    }

    private static boolean implementsMethod(Class<?> aClass, String methodName, Class[] methodParameterTypes) {
        while (true) {
            if (Operator.class.equals(aClass)
                    || !Operator.class.isAssignableFrom(aClass)) {
                return false;
            }
            try {
                Method declaredMethod = aClass.getDeclaredMethod(methodName, methodParameterTypes);
                return declaredMethod.getModifiers() != Modifier.ABSTRACT;
            } catch (NoSuchMethodException e) {
                aClass = aClass.getSuperclass();
            }
        }
    }

    private void initializeOperator() throws OperatorException {
        Assert.state(targetProduct == null, "targetProduct == null");
        Assert.state(operator != null, "operator != null");
        Assert.state(!initialising, "!initialising, attempt to call getTargetProduct() from within initialise()?");

        try {
            initialising = true;
            if (!(operator instanceof GraphOp)) {
                initSourceProductFields();
                injectParameterValues();
                injectConfiguration();
            }
            operator.initialize();
            initTargetProduct();
            initTargetProperties();
            if (!(operator instanceof GraphOp)) {
                initPassThrough();
            }
            initTargetImages();
            initGraphMetadata();

            targetProduct.setModified(false);
        } finally {
            initialising = false;
        }
    }

    /**
     * Updates this operator forcing it to recreate the target product.
     * <i>Warning: Experimental API added by nf (25.02.2010)</i><br/>
     *
     * @since BEAM 4.8
     */
    public void updateOperator() throws OperatorException {
        targetProduct = null;
        initializeOperator();
    }


    private PropertyContainer getOperatorValueContainer() {
        if (propertyContainer == null) {
            PropertyDescriptorFactory parameterDescriptorFactory = new ParameterDescriptorFactory(sourceProductMap);
            propertyContainer = PropertyContainer.createObjectBacked(operator, parameterDescriptorFactory);
        }
        return propertyContainer;
    }

    private void initGraphMetadata() {
        final MetadataElement metadataRoot = targetProduct.getMetadataRoot();
        MetadataElement targetGraphME = metadataRoot.getElement("Processing_Graph");
        if (targetGraphME == null) {
            targetGraphME = new MetadataElement("Processing_Graph");
            metadataRoot.addElement(targetGraphME);
        }
        convertOperatorContextToMetadata(this, targetGraphME);
    }

    private void convertOperatorContextToMetadata(OperatorContext context, MetadataElement targetGraphME) {
        String opId = context.getId();
        boolean contains = false;
        int nodeElementCount = 0;
        for (MetadataElement element : targetGraphME.getElements()) {
            MetadataAttribute idAttribute = element.getAttribute("id");
            if (idAttribute.getData().getElemString().equals(opId)) {
                contains = true;
            }
            if (element.getName().startsWith("node")) {
                nodeElementCount++;
            }
        }
        if (contains) {
            return;
        }
        final String opName = OperatorSpi.getOperatorAlias(context.operator.getClass());
        MetadataElement targetNodeME = new MetadataElement(String.format("node.%d", nodeElementCount));
        targetGraphME.addElement(targetNodeME);
        targetNodeME.addAttribute(new MetadataAttribute("id", ProductData.createInstance(opId), false));
        targetNodeME.addAttribute(new MetadataAttribute("operator", ProductData.createInstance(opName), false));
        final MetadataElement targetSourcesME = new MetadataElement("sources");
        for (String sourceId : context.sourceProductMap.keySet()) {
            final Product sourceProduct = context.sourceProductMap.get(sourceId);
            final String sourceNodeId;
            if (sourceProduct.getFileLocation() != null) {
                sourceNodeId = sourceProduct.getFileLocation().toURI().toASCIIString();
            } else if (sourceProduct.getProductReader() instanceof OperatorProductReader) {
                final OperatorProductReader productReader = (OperatorProductReader) sourceProduct.getProductReader();
                convertOperatorContextToMetadata(productReader.getOperatorContext(), targetGraphME);
                sourceNodeId = productReader.getOperatorContext().getId();
            } else {
                sourceNodeId = "product:" + sourceProduct.getDisplayName();
            }
            final MetadataAttribute sourceAttribute = new MetadataAttribute(sourceId,
                                                                            ProductData.createInstance(sourceNodeId),
                                                                            false);
            targetSourcesME.addAttribute(sourceAttribute);
        }
        targetNodeME.addElement(targetSourcesME);

        final DefaultDomConverter domConverter = new DefaultDomConverter(context.operator.getClass(),
                                                                         new ParameterDescriptorFactory(
                                                                                 sourceProductMap));
        final Xpp3DomElement parametersDom = new Xpp3DomElement("parameters");
        try {
            domConverter.convertValueToDom(context.operator, parametersDom);
        } catch (ConversionException e) {
            e.printStackTrace();
        }
        final MetadataElement targetParametersME = new MetadataElement("parameters");
        addDomToMetadata(parametersDom, targetParametersME);
        targetNodeME.addElement(targetParametersME);
    }

    private static void addDomToMetadata(DomElement parentDE, MetadataElement parentME) {
        final HashMap<String, List<DomElement>> map = new HashMap<String, List<DomElement>>(
                parentDE.getChildCount() + 5);
        for (DomElement childDE : parentDE.getChildren()) {
            final String name = childDE.getName();
            List<DomElement> elementList = map.get(name);
            if (elementList == null) {
                elementList = new ArrayList<DomElement>(3);
                map.put(name, elementList);
            }
            elementList.add(childDE);
        }
        for (Entry<String, List<DomElement>> entry : map.entrySet()) {
            String name = entry.getKey();
            final List<DomElement> elementList = entry.getValue();
            if (elementList.size() > 1) {
                for (int i = 0; i < elementList.size(); i++) {
                    addDomToMetadata(elementList.get(i), name + "." + i, parentME);
                }
            } else {
                addDomToMetadata(elementList.get(0), name, parentME);
            }
        }
    }

    private static void addDomToMetadata(DomElement childDE, String name, MetadataElement parentME) {
        if (childDE.getChildCount() > 0 || childDE.getAttributeNames().length > 0) {
            final MetadataElement childME = new MetadataElement(name);
            addDomToMetadata(childDE, childME);
            parentME.addElement(childME);

            if (childDE.getAttributeNames().length != 0) {
                String[] attributeNames = childDE.getAttributeNames();
                for (String attributeName : attributeNames) {
                    String attributeValue = childDE.getAttribute(attributeName);
                    final ProductData valueMEAtrr = ProductData.createInstance(attributeValue);
                    final MetadataAttribute mdAttribute = new MetadataAttribute(attributeName, valueMEAtrr, true);
                    childME.addAttribute(mdAttribute);
                }
            }
        } else {
            String valueDE = childDE.getValue();
            if (valueDE == null) {
                valueDE = "";
            }
            final ProductData valueME = ProductData.createInstance(valueDE);
            final MetadataAttribute attribute = new MetadataAttribute(name, valueME, true);
            parentME.addAttribute(attribute);
        }
    }

    private boolean initPassThrough() {
        passThrough = false;
        Product[] sourceProducts = getSourceProducts();
        for (Product sourceProduct : sourceProducts) {
            if (targetProduct == sourceProduct) {
                passThrough = true;
            }
        }
        return passThrough;
    }

    private void initTargetImages() {
        if (targetProduct.getPreferredTileSize() == null) {
            targetProduct.setPreferredTileSize(getPreferredTileSize());
        }
        final Band[] targetBands = targetProduct.getBands();
        Object[][] locks = null;
        if (operatorMustComputeTileStack()) {
            Dimension tileSize = targetProduct.getPreferredTileSize();
            int width = targetProduct.getSceneRasterWidth();
            int height = targetProduct.getSceneRasterHeight();
            locks = OperatorImageTileStack.createLocks(width, height, tileSize);
        }
        targetImageMap = new HashMap<Band, OperatorImage>(targetBands.length * 2);
        for (final Band targetBand : targetBands) {
            // Only register non-virtual bands
            if (isRegularBand(targetBand)) {

                OperatorImage opImage = getOwnedOperatorImage(targetBand);
                if (opImage == null) {

                    final OperatorImage image;
                    if (operatorMustComputeTileStack()) {
                        image = new OperatorImageTileStack(targetBand, this, locks);
                    } else {
                        image = new OperatorImage(targetBand, this);
                    }
                    targetImageMap.put(targetBand, image);

                    // Note: It is legal not to set the newly created operator image
                    // in the target band, if it already has a source image set.
                    // This case occurs for "pass-through" operators.
                    // Pull processing in GPF is primarily triggered by fetching tiles
                    // of the target images of this operator context, not directly
                    // by using the band's source images. Otherwise the WriteOp.computeTile()
                    // method would never be called.
                    //
                    if (!targetBand.isSourceImageSet()) {
                        targetBand.setSourceImage(image);
                    }
                } else {
                    targetBand.getSourceImage().reset();
                    targetImageMap.put(targetBand, opImage);
                }
            }

        }
    }

    private OperatorImage getOwnedOperatorImage(Band targetBand) {
        // If there is no source image then we can't own it
        if (!targetBand.isSourceImageSet()) {
            return null;
        }
        // If the source image is not an OperatorImage then we can't own it neither
        RenderedImage renderedImage = targetBand.getSourceImage().getImage(0);
        if (!(renderedImage instanceof OperatorImage)) {
            return null;
        }
        // If the OperatorImage's context is not us, then it is not ours
        OperatorImage operatorImage = (OperatorImage) renderedImage;
        if (this != operatorImage.getOperatorContext()) {
            return null;
        }
        // Now it must be an OperatorImage that we have created
        return operatorImage;
    }

    private Dimension getPreferredTileSize() {
        Dimension tileSize = null;
        for (final Product sourceProduct : sourceProductList) {
            if (sourceProduct.getPreferredTileSize() != null &&
                    sourceProduct.getSceneRasterWidth() == targetProduct.getSceneRasterWidth() &&
                    sourceProduct.getSceneRasterHeight() == targetProduct.getSceneRasterHeight()) {
                tileSize = sourceProduct.getPreferredTileSize();
                break;
            }
        }
        if (tileSize == null) {
            tileSize = JAIUtils.computePreferredTileSize(targetProduct.getSceneRasterWidth(),
                                                         targetProduct.getSceneRasterHeight(), 4);
        }
        return tileSize;
    }

    public static boolean isRegularBand(Band targetBand) {
        // Note: "instanceof" has intentionally not been used here.
        return targetBand.getClass() == Band.class;
    }

    private void initTargetProduct() throws OperatorException {
        Class<? extends Operator> operatorClass = operator.getClass();
        initTargetProduct(operatorClass);
        if (targetProduct == null) {
            final String message = formatExceptionMessage("No target product set.");
            throw new OperatorException(message);
        }
        if (targetProduct.getProductReader() == null) {
            targetProduct.setProductReader(new OperatorProductReader(this));
        }
        if (renderingHints != null && GPF.KEY_TILE_SIZE.isCompatibleValue(renderingHints.get(GPF.KEY_TILE_SIZE))) {
            targetProduct.setPreferredTileSize((Dimension) renderingHints.get(GPF.KEY_TILE_SIZE));
        }
    }

    private void initTargetProduct(Class<? extends Operator> operatorClass) {
        Class<?> superClass = operatorClass.getSuperclass();
        if (superClass != null && !superClass.equals(Operator.class)) {
            initTargetProduct((Class<? extends Operator>) superClass);
        }
        for (Field declaredField : operatorClass.getDeclaredFields()) {
            TargetProduct targetProductAnnotation = declaredField.getAnnotation(TargetProduct.class);
            if (targetProductAnnotation != null) {
                if (!declaredField.getType().equals(Product.class)) {
                    String msg = formatExceptionMessage("Field '%s' annotated as target product is not of type '%s'.",
                                                        declaredField.getName(), Product.class);
                    throw new OperatorException(msg);
                }
                final Product targetProduct = (Product) getOperatorFieldValue(declaredField);
                if (targetProduct != null) {
                    this.targetProduct = targetProduct;
                } else {
                    if (this.targetProduct != null) {
                        setOperatorFieldValue(declaredField, this.targetProduct);
                    } else {
                        final String message = formatExceptionMessage("No target product set.");
                        throw new OperatorException(message);
                    }
                }
            }
        }
    }

    private void initTargetProperties() throws OperatorException {
        Field[] declaredFields = operator.getClass().getDeclaredFields();
        for (Field declaredField : declaredFields) {
            TargetProperty targetPropertyAnnotation = declaredField.getAnnotation(TargetProperty.class);
            if (targetPropertyAnnotation != null) {
                Object propertyValue = getOperatorFieldValue(declaredField);
                String fieldName = declaredField.getName();
                if (targetPropertyMap.containsKey(fieldName)) {
                    final String message = formatExceptionMessage(
                            "Name of field '%s' is already used as target property alias.",
                            fieldName);
                    throw new OperatorException(message);
                }
                targetPropertyMap.put(fieldName, propertyValue);
                if (!targetPropertyAnnotation.alias().isEmpty()) {
                    String aliasName = targetPropertyAnnotation.alias();
                    if (targetPropertyMap.containsKey(aliasName)) {
                        final String message = formatExceptionMessage(
                                "Alias of field '%s' is already used by another target property.",
                                aliasName);
                        throw new OperatorException(message);
                    }
                    targetPropertyMap.put(aliasName, propertyValue);
                }
            }
        }
    }

    private void initSourceProductFields() throws OperatorException {
        initSourceProductFields(operator.getClass());
    }

    private void initSourceProductFields(Class<? extends Operator> operatorClass) {
        Field[] declaredFields = operatorClass.getDeclaredFields();
        for (Field declaredField : declaredFields) {
            SourceProduct sourceProductAnnotation = declaredField.getAnnotation(SourceProduct.class);
            if (sourceProductAnnotation != null) {
                processSourceProductField(declaredField, sourceProductAnnotation);
            }
            SourceProducts sourceProductsAnnotation = declaredField.getAnnotation(SourceProducts.class);
            if (sourceProductsAnnotation != null) {
                processSourceProductsField(declaredField, sourceProductsAnnotation);
            }
        }
        Class<?> superClass = operatorClass.getSuperclass();
        if (superClass != null && !superClass.equals(Operator.class)) {
            initSourceProductFields((Class<? extends Operator>) superClass);
        }
    }

    private void processSourceProductField(Field declaredField, SourceProduct sourceProductAnnotation) throws OperatorException {
        if (declaredField.getType().equals(Product.class)) {
            Product sourceProduct = getSourceProduct(declaredField.getName());
            if (sourceProduct == null) {
                sourceProduct = getSourceProduct(sourceProductAnnotation.alias());
            }
            if (sourceProduct != null) {
                validateSourceProduct(declaredField.getName(),
                                      sourceProduct,
                                      sourceProductAnnotation.type(),
                                      sourceProductAnnotation.bands());
                setSourceProductFieldValue(declaredField, sourceProduct);
            } else {
                sourceProduct = getSourceProductFieldValue(declaredField);
                if (sourceProduct != null) {
                    setSourceProduct(declaredField.getName(), sourceProduct);
                } else if (!sourceProductAnnotation.optional()) {
                    String text = "Mandatory source product (field '%s') not set.";
                    String msg = formatExceptionMessage(text, declaredField.getName());
                    throw new OperatorException(msg);
                }
            }
        } else {
            String text = "A source product (field '%s') must be of type '%s'.";
            String msg = formatExceptionMessage(text, declaredField.getName(), Product.class.getName());
            throw new OperatorException(msg);
        }
    }

    private void processSourceProductsField(Field declaredField, SourceProducts sourceProductsAnnotation) throws
            OperatorException {
        if (declaredField.getType().equals(Product[].class)) {
            Product[] sourceProducts = getSourceProducts();
            if (sourceProducts.length > 0) {
                setSourceProductsFieldValue(declaredField, sourceProducts);//getUnnamedProducts());
            } else {
                sourceProducts = getSourceProductsFieldValue(declaredField);
                if (sourceProducts != null) {
                    for (int i = 0; i < sourceProducts.length; i++) {
                        Product sourceProduct = sourceProducts[i];
                        setSourceProduct(GPF.SOURCE_PRODUCT_FIELD_NAME + i, sourceProduct);
                    }
                }
                sourceProducts = getSourceProducts();
            }
            if (sourceProductsAnnotation.count() < 0) {
                if (sourceProducts.length == 0) {
                    String msg = formatExceptionMessage("At least a single source product expected.");
                    throw new OperatorException(msg);
                }
            } else if (sourceProductsAnnotation.count() > 0) {
                if (sourceProductsAnnotation.count() != sourceProducts.length) {
                    String text = "Wrong number of source products. Required %d, found %d.";
                    String msg = formatExceptionMessage(text, sourceProductsAnnotation.count(), sourceProducts.length);
                    throw new OperatorException(msg);
                }
            }
            for (Product sourceProduct : sourceProducts) {
                validateSourceProduct(declaredField.getName(),
                                      sourceProduct,
                                      sourceProductsAnnotation.type(),
                                      sourceProductsAnnotation.bands());
            }
        } else {
            String text = "Source products (field '%s') must be of type '%s'.";
            String msg = formatExceptionMessage(text, declaredField.getName(), Product[].class.getName());
            throw new OperatorException(msg);
        }
    }

    private Product[] getUnnamedProducts() {
        final Map<String, Product> map = new HashMap<String, Product>(sourceProductMap);
        final Field[] sourceProductFields = getAnnotatedSourceProductFields(operator);
        for (Field sourceProductField : sourceProductFields) {
            final SourceProduct annotation = sourceProductField.getAnnotation(SourceProduct.class);
            map.remove(sourceProductField.getName());
            map.remove(annotation.alias());
        }
        final Collection<Product> unnamedProductList = map.values();
        return unnamedProductList.toArray(new Product[unnamedProductList.size()]);
    }

    private Field[] getAnnotatedSourceProductFields(Operator operator1) {
        Field[] declaredFields = operator1.getClass().getDeclaredFields();
        List<Field> fieldList = new ArrayList<Field>();
        for (Field declaredField : declaredFields) {
            SourceProduct sourceProductAnnotation = declaredField.getAnnotation(SourceProduct.class);
            if (sourceProductAnnotation != null) {
                fieldList.add(declaredField);
            }
        }
        return fieldList.toArray(new Field[fieldList.size()]);
    }

    private Product getSourceProductFieldValue(Field declaredField) throws OperatorException {
        return (Product) getOperatorFieldValue(declaredField);
    }

    private void setSourceProductFieldValue(Field declaredField, Product sourceProduct) throws OperatorException {
        setOperatorFieldValue(declaredField, sourceProduct);
    }

    private Product[] getSourceProductsFieldValue(Field declaredField) throws OperatorException {
        return (Product[]) getOperatorFieldValue(declaredField);
    }

    private void setSourceProductsFieldValue(Field declaredField, Product[] sourceProducts) throws OperatorException {
        setOperatorFieldValue(declaredField, sourceProducts);
    }

    private void validateSourceProduct(String fieldName, Product sourceProduct, String typeRegExp,
                                       String[] bandNames) throws OperatorException {
        if (!typeRegExp.isEmpty()) {
            final String productType = sourceProduct.getProductType();
            if (!typeRegExp.equalsIgnoreCase(productType) && !Pattern.matches(typeRegExp, productType)) {
                String msg = formatExceptionMessage(
                        "A source product (field '%s') of type '%s' does not match type '%s'",
                        fieldName,
                        productType,
                        typeRegExp);
                throw new OperatorException(msg);
            }
        }

        if (bandNames.length != 0) {
            for (String bandName : bandNames) {
                if (!sourceProduct.containsBand(bandName)) {
                    String msg = formatExceptionMessage("A source product (field '%s') does not contain the band '%s'",
                                                        fieldName, bandName);
                    throw new OperatorException(msg);
                }
            }
        }
    }


    private Object getOperatorFieldValue(Field declaredField) throws OperatorException {
        boolean oldState = declaredField.isAccessible();
        try {
            declaredField.setAccessible(true);
            try {
                return declaredField.get(getOperator());
            } catch (IllegalAccessException e) {
                String msg = formatExceptionMessage("Unable to get declared field '%s'.", declaredField.getName());
                throw new OperatorException(msg, e);
            }
        } finally {
            declaredField.setAccessible(oldState);
        }
    }

    private void setOperatorFieldValue(Field declaredField, Object value) throws OperatorException {
        boolean oldState = declaredField.isAccessible();
        try {
            declaredField.setAccessible(true);
            try {
                declaredField.set(getOperator(), value);
            } catch (IllegalAccessException e) {
                String msg = formatExceptionMessage("Unable to set declared field '%s'", declaredField.getName());
                throw new OperatorException(msg, e);
            }
        } finally {
            declaredField.setAccessible(oldState);
        }
    }

    public void injectConfiguration() throws OperatorException {
        if (configuration != null) {
            try {
                configureOperator(operator, configuration);
            } catch (OperatorException e) {
                throw e;
            } catch (Throwable t) {
                throw new OperatorException(formatExceptionMessage("%s", t.getMessage()), t);
            }
        }
    }

    private void configureOperator(Operator operator, OperatorConfiguration operatorConfiguration) throws
            ValidationException,
            ConversionException {
        ParameterDescriptorFactory parameterDescriptorFactory = new ParameterDescriptorFactory(sourceProductMap);
        DefaultDomConverter domConverter = new DefaultDomConverter(operator.getClass(), parameterDescriptorFactory);
        domConverter.convertDomToValue(operatorConfiguration.getConfiguration(), operator);
        PropertyContainer propertyContainer = PropertyContainer.createObjectBacked(operator, parameterDescriptorFactory);
        Set<Reference> referenceSet = operatorConfiguration.getReferenceSet();
        for (Reference reference : referenceSet) {
            Property property = propertyContainer.getProperty(reference.getParameterName());
            property.setValue(reference.getValue());
        }
    }

    public void injectParameterDefaultValues() throws OperatorException {
        try {
            getOperatorValueContainer().setDefaultValues();
        } catch (ValidationException e) {
            throw new OperatorException(formatExceptionMessage("%s", e.getMessage()), e);
        }
    }

    private void injectParameterValues() throws OperatorException {
        if (parameters != null) {
            for (String parameterName : parameters.keySet()) {
                final Property property = getOperatorValueContainer().getProperty(parameterName);
                if (property == null) {
                    throw new OperatorException(formatExceptionMessage("Unknown parameter '%s'.", parameterName));
                }
                try {
                    PropertyDescriptor descriptor = property.getDescriptor();
                    if (descriptor.getAttribute(RasterDataNodeValues.ATTRIBUTE_NAME) != null) {
                        Product sourceProduct = sourceProductList.get(0);
                        if (sourceProduct == null) {
                            throw new OperatorException(formatExceptionMessage("No source produt."));
                        }
                        Object object = descriptor.getAttribute(RasterDataNodeValues.ATTRIBUTE_NAME);
                        Class<? extends RasterDataNode> rasterDataNodeType = (Class<? extends RasterDataNode>) object;
                        String[] names = RasterDataNodeValues.getNames(sourceProduct, rasterDataNodeType);
                        ValueSet valueSet = new ValueSet(names);
                        descriptor.setValueSet(valueSet);
                    }
                    property.setValue(parameters.get(parameterName));
                } catch (ValidationException e) {
                    throw new OperatorException(formatExceptionMessage("%s", e.getMessage()), e);
                }
            }
        }
    }

    private String formatExceptionMessage(String format, Object... args) {
        Object[] allArgs = new Object[args.length + 1];
        allArgs[0] = operator.getClass().getSimpleName();
        System.arraycopy(args, 0, allArgs, 1, allArgs.length - 1);
        return String.format("Operator '%s': " + format, allArgs);
    }

    /**
     * Sums the target nanos/pixel of all source images.
     *
     * @return the nanos/pixel spend in the sources of this operator.
     */
    double getSourceNanosPerPixel() {
        double sourceNanosPerPixel = 0;
        for (Product product : sourceProductList) {
            OperatorContext operatorContext = getOperatorContext(product);
            if (operatorContext != null) {
                sourceNanosPerPixel += operatorContext.getPerformanceMetric().getTargetNanosPerPixel();
            }
        }
        return sourceNanosPerPixel;
    }

    /**
     * @param product The product.
     * @return The operator context for the given product, or null.
     */
    static OperatorContext getOperatorContext(Product product) {
        ProductReader productReader = product.getProductReader();
        if (productReader instanceof OperatorProductReader) {
            OperatorProductReader operatorProductReader = (OperatorProductReader) productReader;
            return operatorProductReader.getOperatorContext();
        }
        return null;
    }

    public void logPerformanceAnalysis() {
        logPerformanceAnalysis("");
    }

    private void logPerformanceAnalysis(String indent) {
        logger.info(indent + "Performance analysis for operator " + getOperatorSpi().getOperatorAlias() + ":");
        logger.info(indent + "  Net: " + getPerformanceMetric().getNanosPerPixel() + " ns/pixel");
        logger.info(indent + "  Target: " + getPerformanceMetric().getTargetNanosPerPixel() + " ns/pixel");
        logger.info(
                indent + "  Sources (" + sourceProductList.size() + "): " + getPerformanceMetric().getSourceNanosPerPixel() + " ns/pixel");
        for (Product product : sourceProductList) {
            OperatorContext operatorContext = getOperatorContext(product);
            if (operatorContext != null) {
                operatorContext.logPerformanceAnalysis(indent + "  ");
            }
        }
    }

    boolean isComputing(Band band) {
        if (band.isSourceImageSet()) {
            RenderedImage sourceImage = band.getSourceImage().getImage(0);
            OperatorImage targetImage = getTargetImage(band);
            return targetImage == sourceImage;
        } else {
            return false;
        }
    }

    public boolean requiresAllBands() {
        return requiresAllBands;
    }

    public void setRequiresAllBands(boolean requiresAllBands) {
        this.requiresAllBands = requiresAllBands;
    }
}
