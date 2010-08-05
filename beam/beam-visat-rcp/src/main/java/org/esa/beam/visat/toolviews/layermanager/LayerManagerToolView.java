/*
 * Copyright (C) 2010 Brockmann Consult GmbH (info@brockmann-consult.de)
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
package org.esa.beam.visat.toolviews.layermanager;

import com.bc.ceres.glayer.Layer;
import org.esa.beam.framework.ui.product.ProductSceneView;

import javax.swing.JComponent;
import javax.swing.JPanel;
import java.awt.BorderLayout;
import java.util.WeakHashMap;

/**
 * Layer manager tool view.
 * <p/>
 * <i>Note: This API is not public yet and may significantly change in the future. Use it at your own risk.</i>
 *
 * @author Norman Fomferra
 * @version $Revision: 1.5 $ $Date: 2010-08-05 17:00:56 $
 * @since BEAM 4.6
 */
public class LayerManagerToolView extends AbstractLayerToolView {

    private WeakHashMap<ProductSceneView, LayerManagerForm> layerManagerMap;
    private LayerManagerForm activeForm;

    public LayerManagerToolView() {
    }

    @Override
    protected JComponent createControl() {
        layerManagerMap = new WeakHashMap<ProductSceneView, LayerManagerForm>();
        return super.createControl();
    }

    @Override
    protected void viewClosed(ProductSceneView view) {
        layerManagerMap.remove(view);
    }

    @Override
    protected void viewSelectionChanged(ProductSceneView oldView, ProductSceneView newView) {
        realizeActiveForm();
    }

    @Override
    protected void layerSelectionChanged(Layer oldLayer, Layer selectedLayer) {
        if (activeForm != null) {
            activeForm.updateFormControl();
        }
    }

    private void realizeActiveForm() {
        final JPanel controlPanel = getControlPanel();

        if (controlPanel.getComponentCount() > 0) {
            controlPanel.remove(0);
        }

        if (getSelectedView() != null) {
            activeForm = getOrCreateActiveForm(getSelectedView());
            controlPanel.add(activeForm.getFormControl(), BorderLayout.CENTER);
        } else {
            activeForm = null;
        }

        controlPanel.validate();
        controlPanel.repaint();
    }

    protected LayerManagerForm getOrCreateActiveForm(ProductSceneView view) {
        if (layerManagerMap.containsKey(view)) {
            activeForm = layerManagerMap.get(view);
        } else {
            activeForm = new LayerManagerForm(getAppContext(), getDescriptor().getHelpId());
            layerManagerMap.put(view, activeForm);
        }
        return activeForm;
    }

}
