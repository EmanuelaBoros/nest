/*
 * $Id: CopyMaskDialog.java,v 1.1 2009-12-04 19:06:45 lveci Exp $
 *
 * Copyright (C) 2009 by Brockmann Consult (info@brockmann-consult.de)
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
package org.esa.beam.visat.toolviews.mask;

import com.bc.ceres.swing.TableLayout;

import org.esa.beam.framework.datamodel.GeoCoding;
import org.esa.beam.framework.datamodel.Mask;
import org.esa.beam.framework.datamodel.Product;
import org.esa.beam.framework.ui.ModalDialog;
import org.esa.beam.util.ProductUtils;

import java.awt.Rectangle;
import java.awt.geom.GeneralPath;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.swing.ButtonGroup;
import javax.swing.ButtonModel;
import javax.swing.JCheckBox;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.SwingConstants;

/**
 * @author Marco Zuehlke
 * @version $Revision: 1.1 $ $Date: 2009-12-04 19:06:45 $
 * @since BEAM 4.7
 */
class CopyMaskDialog extends ModalDialog {

    private final Mask[] selectedMasks;
    private final Product sourceProduct;
    private final Product[] allProducts;
    private final Map<Product, ButtonModel> definitionMap;
    private final Map<Product, ButtonModel> dataMap;

    CopyMaskDialog(Product product, Product[] allProducts, Mask[] selectedMasks) {
        super(null, "Copy Mask(s)", ModalDialog.ID_OK_CANCEL | ModalDialog.ID_HELP, "copyMaskEditor");
        this.sourceProduct = product;
        this.allProducts = allProducts;
        this.selectedMasks = selectedMasks;
        definitionMap = new HashMap<Product, ButtonModel>();
        dataMap = new HashMap<Product, ButtonModel>();
        getJDialog().setResizable(false);
        setContent(createUI());
    }
    
    Product[] getMaskPixelTargets() {
        return getSelectedProducts(dataMap);
    }
    
    Product[] getMaskDefinitionTargets() {
        return getSelectedProducts(definitionMap);
    }
    
    private static Product[] getSelectedProducts(Map<Product, ButtonModel> buttonMap) {
        List<Product> selectedProducts = new ArrayList<Product>(buttonMap.size());
        for (Map.Entry<Product, ButtonModel> entry : buttonMap.entrySet()) {
            Product product = entry.getKey();
            ButtonModel buttonModel = entry.getValue();
            if (buttonModel.isSelected()) {
                selectedProducts.add(product);
            }
        }
        return selectedProducts.toArray(new Product[selectedProducts.size()]);
    }
    
    private JComponent createUI() {
        final TableLayout layout = new TableLayout(3);
        layout.setTableFill(TableLayout.Fill.BOTH);
        layout.setTableAnchor(TableLayout.Anchor.WEST);
        layout.setTableWeightX(1.0);
        layout.setTableWeightY(0.0);
        layout.setTablePadding(5, 5);
        
        
        final JPanel panel = new JPanel(layout);
        panel.add(new JLabel(" "));
        panel.add(new JLabel("Definition"));
        panel.add(new JLabel("Pixels"));
        int row = 1;
        for (Product targetProduct : allProducts) {
            if (targetProduct != sourceProduct) {
                panel.add(new JLabel(targetProduct.getDisplayName()));
                
                boolean canCopyDef = canCopyDefinition(selectedMasks, targetProduct);
                JCheckBox defCheckBox = createCeckbox(panel, canCopyDef);
                if (canCopyDef) {
                    definitionMap.put(targetProduct, defCheckBox.getModel());
                }
                
                boolean canCopyData = intersectsWith(sourceProduct, targetProduct);
                JCheckBox dataCheckBox = createCeckbox(panel, canCopyData);
                if (canCopyData) {
                    dataMap.put(targetProduct, dataCheckBox.getModel());
                }
                
                if (canCopyData && canCopyDef) {
                    ButtonGroup buttonGroup = new ButtonGroup();
                    buttonGroup.add(dataCheckBox);
                    buttonGroup.add(defCheckBox);
                }
                row++;
            }
        }
        return panel;
    }

    private static JCheckBox createCeckbox(final JPanel panel, boolean enabled) {
        JCheckBox checkBox = new JCheckBox();
        checkBox.setHorizontalAlignment(SwingConstants.CENTER);
        checkBox.setEnabled(enabled);
        panel.add(checkBox);
        return checkBox;
    }

    private static boolean canCopyDefinition(Mask[] masks, Product targetProduct) {
        boolean canCopyDef = true;
        for (Mask mask : masks) {
            boolean canTransferMask = mask.getImageType().canTransferMask(mask, targetProduct);
            canCopyDef = canCopyDef && canTransferMask;
        }
        return canCopyDef;
    }
    
    private static boolean intersectsWith(Product sourceProduct, Product targetProduct) {
        final GeoCoding srcGeoCoding = sourceProduct.getGeoCoding();
        final GeoCoding targetGeoCoding = targetProduct.getGeoCoding();
        if (srcGeoCoding.canGetGeoPos() && targetGeoCoding.canGetGeoPos()) {
            final GeneralPath[] sourcePath = ProductUtils.createGeoBoundaryPaths(sourceProduct);
            final GeneralPath[] targetPath = ProductUtils.createGeoBoundaryPaths(targetProduct);
            for (GeneralPath spath : sourcePath) {
                Rectangle bounds = spath.getBounds();
                for (GeneralPath tPath : targetPath) {
                    if (tPath.getBounds().intersects(bounds)) {
                        return true;
                    }
                }
            }
        }
        return false;
    }
}
