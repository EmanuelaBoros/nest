package org.esa.beam.visat.toolviews.layermanager;

import org.esa.beam.framework.ui.AppContext;
import org.esa.beam.framework.ui.UIUtils;

import javax.swing.AbstractAction;
import javax.swing.Action;
import java.awt.event.ActionEvent;

import com.bc.ceres.glayer.Layer;

/**
*
* @author Marco Peters
* @version $Revision: 1.1 $ $Date: 2009-04-27 13:08:25 $
* @since BEAM 4.6
*/
class MoveLayerDownAction extends AbstractAction {

    private final AppContext appContext;

    MoveLayerDownAction(AppContext appContext) {
        super("Move Layer Down", UIUtils.loadImageIcon("icons/Down24.gif"));
        this.appContext = appContext;
        putValue(Action.ACTION_COMMAND_KEY, getClass().getName());
    }


    @Override
    public void actionPerformed(ActionEvent e) {
        final Layer selectedLayer = appContext.getSelectedProductSceneView().getSelectedLayer();
        Layer rootLayer = appContext.getSelectedProductSceneView().getRootLayer();
        if (selectedLayer != null && rootLayer != selectedLayer) {
            moveDown(selectedLayer);
        }
    }

    void moveDown(Layer layer) {
        final Layer parentLayer = layer.getParent();
        final int layerIndex = parentLayer.getChildIndex(layer.getId());

        final boolean isLast = layerIndex == (parentLayer.getChildren().size() - 1);
        if (!isLast) {
            parentLayer.getChildren().remove(layer);
            parentLayer.getChildren().add(layerIndex + 1, layer);
        }
    }
}