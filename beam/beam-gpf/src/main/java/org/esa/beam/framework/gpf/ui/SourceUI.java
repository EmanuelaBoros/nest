package org.esa.beam.framework.gpf.ui;

import org.esa.beam.framework.ui.AppContext;
import org.esa.beam.framework.datamodel.Product;
import org.esa.beam.framework.dataio.ProductIO;

import javax.swing.*;
import java.util.ArrayList;
import java.util.Map;
import java.util.List;
import java.io.IOException;
import java.io.File;

import com.bc.ceres.swing.TableLayout;


/**
 * Reader OperatorUI
 */
public class SourceUI extends BaseOperatorUI {

    SourceProductSelector sourceProductSelector = null;
    private static final String FILE_PARAMETER = "file";

    @Override
    public JComponent CreateOpTab(String operatorName, Map<String, Object> parameterMap, AppContext appContext) {

        paramMap = parameterMap;
        final List<SourceProductSelector> sourceProductSelectorList = new ArrayList<SourceProductSelector>(3);
        sourceProductSelector = new SourceProductSelector(appContext);
        sourceProductSelectorList.add(sourceProductSelector);

        final TableLayout tableLayout = new TableLayout(1);
        tableLayout.setTableAnchor(TableLayout.Anchor.WEST);
        tableLayout.setTableWeightX(1.0);
        tableLayout.setTableFill(TableLayout.Fill.HORIZONTAL);
        tableLayout.setTablePadding(3, 3);

        final JPanel ioParametersPanel = new JPanel(tableLayout);
        for (SourceProductSelector selector : sourceProductSelectorList) {
            ioParametersPanel.add(selector.createDefaultPanel());
        }
        ioParametersPanel.add(tableLayout.createVerticalSpacer());

        initSourceProductSelectors(sourceProductSelectorList);

        initParameters();
        
        return ioParametersPanel;
    }

     private static void initSourceProductSelectors(java.util.List<SourceProductSelector> sourceProductSelectorList) {
        for (SourceProductSelector sourceProductSelector : sourceProductSelectorList) {
            sourceProductSelector.initProducts();
            if (sourceProductSelector.getProductCount() > 0) {
                sourceProductSelector.setSelectedIndex(0);
            }
        }
    }

    @Override
    public void initParameters() {
        assert(paramMap != null);
        final Object value = paramMap.get(FILE_PARAMETER);
        if(value != null) {

            try {
                final Product product = ProductIO.readProduct((File)value, null);
                sourceProductSelector.setSelectedProduct(product);
            } catch (IOException e) {
                // do nothing
            }
        }
    }

    @Override
    public UIValidation validateParameters() {
        if(sourceProductSelector != null) {
            if(sourceProductSelector.getSelectedProduct() == null)
                return new UIValidation(UIValidation.State.ERROR, "Source product not selected");
        }
        return new UIValidation(UIValidation.State.OK, "");
    }

    @Override
    public void updateParameters() {
        if(sourceProductSelector != null) {
            final Product prod = sourceProductSelector.getSelectedProduct();
            if(prod != null && prod.getFileLocation() != null)
                paramMap.put(FILE_PARAMETER, prod.getFileLocation());
        }
    }
}
