/*
 * Copyright (C) 2011 Brockmann Consult GmbH (info@brockmann-consult.de)
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

package org.esa.beam.framework.gpf.ui;

import com.bc.ceres.binding.Property;
import com.bc.ceres.binding.PropertyDescriptor;
import com.bc.ceres.binding.PropertySet;
import com.bc.ceres.binding.ValueSet;
import com.bc.ceres.swing.binding.PropertyPane;
import com.bc.ceres.swing.selection.AbstractSelectionChangeListener;
import com.bc.ceres.swing.selection.Selection;
import com.bc.ceres.swing.selection.SelectionChangeEvent;
import org.esa.beam.framework.datamodel.Product;
import org.esa.beam.framework.datamodel.ProductNodeEvent;
import org.esa.beam.framework.datamodel.ProductNodeListener;
import org.esa.beam.framework.datamodel.RasterDataNode;
import org.esa.beam.framework.gpf.GPF;
import org.esa.beam.framework.gpf.OperatorSpi;
import org.esa.beam.framework.gpf.internal.RasterDataNodeValues;
import org.esa.beam.framework.ui.AppContext;

import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTabbedPane;
import javax.swing.border.EmptyBorder;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

// todo (mp, 2008/04/22) add abillity to set the ProductFilter to SourceProductSelectors

/**
 * WARNING: This class belongs to a preliminary API and may change in future releases.
 *
 * @author Norman Fomferra
 * @version $Revision: 1.18 $ $Date: 2011-06-16 19:12:39 $
 */
public class DefaultSingleTargetProductDialog extends SingleTargetProductDialog {

    protected final String operatorName;
    protected final JTabbedPane form;
    protected PropertyDescriptor[] rasterDataNodeTypeProperties;
    protected String targetProductNameSuffix;
    protected ProductChangedHandler productChangedHandler;
    protected DefaultIOParametersPanel ioParametersPanel;
    protected final OperatorParameterSupport parameterSupport;

    public static SingleTargetProductDialog createDefaultDialog(String operatorName, AppContext appContext) {
        return new DefaultSingleTargetProductDialog(operatorName, appContext, operatorName, null);
    }

    public DefaultSingleTargetProductDialog(String operatorName, AppContext appContext, String title, String helpID) {
        super(appContext, title, ID_APPLY_CLOSE_HELP, helpID);
        this.operatorName = operatorName;
        targetProductNameSuffix = "";

        final OperatorSpi operatorSpi = GPF.getDefaultInstance().getOperatorSpiRegistry().getOperatorSpi(operatorName);
        if (operatorSpi == null) {
            throw new IllegalArgumentException("operatorName");
        }

        ioParametersPanel = new DefaultIOParametersPanel(getAppContext(), operatorSpi, getTargetProductSelector());

        this.form = new JTabbedPane();
        this.form.add("I/O Parameters", ioParametersPanel);

        parameterSupport = new OperatorParameterSupport(operatorSpi.getOperatorClass());

        addParametersPane(operatorSpi, appContext, helpID);

        final ArrayList<SourceProductSelector> sourceProductSelectorList = ioParametersPanel.getSourceProductSelectorList();
        if (!sourceProductSelectorList.isEmpty()) {
            productChangedHandler = new ProductChangedHandler();
            sourceProductSelectorList.get(0).addSelectionChangeListener(productChangedHandler);
        }
    }

    protected void addParametersPane(final OperatorSpi operatorSpi, final AppContext appContext, final String helpID) {
        //OperatorMenu operatorMenu = new OperatorMenu(this.getJDialog(),
        //                                             operatorSpi.getOperatorClass(),
        //                                             parameterSupport,
        //                                             helpID);
        PropertySet propertyContainer = parameterSupport.getPopertySet();
        final ArrayList<SourceProductSelector> sourceProductSelectorList = ioParametersPanel.getSourceProductSelectorList();
        if (propertyContainer.getProperties().length > 0) {
            if (!sourceProductSelectorList.isEmpty()) {
                Property[] properties = propertyContainer.getProperties();
                List<PropertyDescriptor> rdnTypeProperties = new ArrayList<PropertyDescriptor>(properties.length);
                for (Property property : properties) {
                    PropertyDescriptor parameterDescriptor = property.getDescriptor();
                    if (parameterDescriptor.getAttribute(RasterDataNodeValues.ATTRIBUTE_NAME) != null) {
                        rdnTypeProperties.add(parameterDescriptor);
                    }
                }
                rasterDataNodeTypeProperties = rdnTypeProperties.toArray(
                        new PropertyDescriptor[rdnTypeProperties.size()]);
            }
            PropertyPane parametersPane = new PropertyPane(propertyContainer);
            final JPanel parametersPanel = parametersPane.createPanel();
            parametersPanel.setBorder(new EmptyBorder(4, 4, 4, 4));
            this.form.add("Processing Parameters", new JScrollPane(parametersPanel));

            //getJDialog().setJMenuBar(operatorMenu.createDefaultMenu());
        }
    }

    @Override
    public int show() {
        ioParametersPanel.initSourceProductSelectors();
        setContent(form);
        return super.show();
    }

    @Override
    public void hide() {
        productChangedHandler.releaseProduct();
        ioParametersPanel.releaseSourceProductSelectors();
        super.hide();
    }

    @Override
    protected Product createTargetProduct() throws Exception {
        final HashMap<String, Product> sourceProducts = ioParametersPanel.createSourceProductsMap();
        return GPF.createProduct(operatorName, parameterSupport.getParameterMap(), sourceProducts);
    }

    public String getTargetProductNameSuffix() {
        return targetProductNameSuffix;
    }

    public void setTargetProductNameSuffix(String suffix) {
        targetProductNameSuffix = suffix;
    }

    private class ProductChangedHandler extends AbstractSelectionChangeListener implements ProductNodeListener {

        private Product currentProduct;

        public void releaseProduct() {
            if (currentProduct != null) {
                currentProduct.removeProductNodeListener(this);
                currentProduct = null;
            }
        }

        @Override
        public void selectionChanged(SelectionChangeEvent event) {
            Selection selection = event.getSelection();
            if (selection != null) {
                final Product selectedProduct = (Product) selection.getSelectedValue();
                if (selectedProduct != currentProduct) {
                    if (currentProduct != null) {
                        currentProduct.removeProductNodeListener(this);
                    }
                    currentProduct = selectedProduct;
                    if (currentProduct != null) {
                        currentProduct.addProductNodeListener(this);
                    }
                    updateTargetProductname();
                    updateValueSets(currentProduct);
                }
            }
        }

        @Override
        public void nodeAdded(ProductNodeEvent event) {
            handleProductNodeEvent(event);
        }

        @Override
        public void nodeChanged(ProductNodeEvent event) {
            handleProductNodeEvent(event);
        }

        @Override
        public void nodeDataChanged(ProductNodeEvent event) {
            handleProductNodeEvent(event);
        }

        @Override
        public void nodeRemoved(ProductNodeEvent event) {
            handleProductNodeEvent(event);
        }

        private void updateTargetProductname() {
            String productName = "";
            if (currentProduct != null) {
                productName = currentProduct.getName();
            }
            final TargetProductSelectorModel targetProductSelectorModel = getTargetProductSelector().getModel();
            targetProductSelectorModel.setProductName(productName + getTargetProductNameSuffix());
        }

        private void handleProductNodeEvent(ProductNodeEvent event) {
            updateValueSets(currentProduct);
        }

        private void updateValueSets(Product product) {
            if (rasterDataNodeTypeProperties != null) {
                for (PropertyDescriptor propertyDescriptor : rasterDataNodeTypeProperties) {
                    updateValueSet(propertyDescriptor, product);
                }
            }
        }
    }

    private static void updateValueSet(PropertyDescriptor propertyDescriptor, Product product) {
        String[] values = new String[0];
        if (product != null) {
            Object object = propertyDescriptor.getAttribute(RasterDataNodeValues.ATTRIBUTE_NAME);
            if (object != null) {
                Class<? extends RasterDataNode> rasterDataNodeType = (Class<? extends RasterDataNode>) object;
                boolean includeEmptyValue = !propertyDescriptor.isNotNull() && !propertyDescriptor.getType().isArray();
                values = RasterDataNodeValues.getNames(product, rasterDataNodeType, includeEmptyValue);
            }
        }
        propertyDescriptor.setValueSet(new ValueSet(values));
    }
}
