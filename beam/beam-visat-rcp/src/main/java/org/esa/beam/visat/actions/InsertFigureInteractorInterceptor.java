package org.esa.beam.visat.actions;

import com.bc.ceres.glayer.Layer;
import com.bc.ceres.glayer.LayerFilter;
import com.bc.ceres.glayer.support.LayerUtils;
import com.bc.ceres.swing.figure.AbstractInteractorListener;
import com.bc.ceres.swing.figure.Interactor;
import org.esa.beam.framework.datamodel.VectorDataNode;
import org.esa.beam.framework.ui.ModalDialog;
import org.esa.beam.framework.ui.product.ProductSceneView;
import org.esa.beam.framework.ui.product.VectorDataLayer;
import org.esa.beam.framework.ui.product.VectorDataLayerFilterFactory;

import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.SwingUtilities;
import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.event.InputEvent;
import java.util.List;


public class InsertFigureInteractorInterceptor extends AbstractInteractorListener {

    @Override
    public boolean canStartInteraction(Interactor interactor, InputEvent inputEvent) {
        ProductSceneView productSceneView = getProductSceneView(inputEvent);
        if (productSceneView == null) {
            return false;
        }

        LayerFilter geometryFilter = VectorDataLayerFilterFactory.createGeometryFilter();

        Layer layer = productSceneView.getSelectedLayer();
        if (geometryFilter.accept(layer)) {
            layer.setVisible(true);
            return true;
        }

        List<Layer> layers = LayerUtils.getChildLayers(productSceneView.getRootLayer(),
                                                       LayerUtils.SEARCH_DEEP, geometryFilter
        );

        VectorDataLayer vectorDataLayer;
        if (layers.isEmpty()) {
            NewVectorDataNodeAction nodeAction = new NewVectorDataNodeAction();
            VectorDataNode vectorDataNode = nodeAction.run();
            LayerFilter nodeFilter = VectorDataLayerFilterFactory.createNodeFilter(vectorDataNode);
            vectorDataLayer = (VectorDataLayer) LayerUtils.getChildLayer(productSceneView.getRootLayer(),
                                                                         LayerUtils.SEARCH_DEEP, nodeFilter);
        } else if (layers.size() == 1) {
            vectorDataLayer = (VectorDataLayer) layers.get(0);
        } else {
            vectorDataLayer = showSelectLayerDialog(productSceneView, layers);
        }
        if (vectorDataLayer == null) {
            // = Cancel
            return false;
        }
        productSceneView.setSelectedLayer(vectorDataLayer);
        if (productSceneView.getSelectedLayer() == vectorDataLayer) {
            vectorDataLayer.setVisible(true);
            return true;
        }
        return false;
    }

    private VectorDataLayer showSelectLayerDialog(ProductSceneView productSceneView, List<Layer> layers) {
        String[] layerNames = new String[layers.size()];
        for (int i = 0; i < layerNames.length; i++) {
            layerNames[i] = layers.get(i).getName();
        }
        JList listBox = new JList(layerNames);
        JPanel panel = new JPanel(new BorderLayout(4, 4));
        panel.add(new JLabel("Please select a geometry container:"), BorderLayout.NORTH);
        panel.add(new JScrollPane(listBox), BorderLayout.CENTER);
        ModalDialog dialog = new ModalDialog(SwingUtilities.getWindowAncestor(productSceneView),
                                             "Select Geometry Container",
                                             ModalDialog.ID_OK_CANCEL_HELP, "");
        dialog.setContent(panel);
        int i = dialog.show();
        if (i == ModalDialog.ID_OK) {
            final int index = listBox.getSelectedIndex();
            if (index >= 0) {
                return (VectorDataLayer) layers.get(index);
            }
        }
        return null;
    }

    private ProductSceneView getProductSceneView(InputEvent event) {
        ProductSceneView productSceneView = null;
        Component component = event.getComponent();
        while (component != null) {
            if (component instanceof ProductSceneView) {
                productSceneView = (ProductSceneView) component;
                break;
            }
            component = component.getParent();
        }
        return productSceneView;
    }

}
