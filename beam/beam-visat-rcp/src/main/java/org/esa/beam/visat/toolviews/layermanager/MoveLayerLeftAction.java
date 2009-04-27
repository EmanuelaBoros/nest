package org.esa.beam.visat.toolviews.layermanager;

import org.esa.beam.framework.ui.AppContext;
import org.esa.beam.framework.ui.UIUtils;

import javax.swing.AbstractAction;
import javax.swing.Action;
import java.awt.event.ActionEvent;

import com.bc.ceres.glayer.Layer;

/**
 * @author Marco Peters
 * @version $Revision: 1.1 $ $Date: 2009-04-27 13:08:25 $
 * @since BEAM 4.6
 */
class MoveLayerLeftAction extends AbstractAction {

    private final AppContext appContext;

    MoveLayerLeftAction(AppContext appContext) {
        super("Move Layer Left", UIUtils.loadImageIcon("icons/Left24.gif"));
        this.appContext = appContext;
        putValue(Action.ACTION_COMMAND_KEY, getClass().getName());
    }


    @Override
    public void actionPerformed(ActionEvent e) {
        final Layer selectedLayer = appContext.getSelectedProductSceneView().getSelectedLayer();
        Layer rootLayer = appContext.getSelectedProductSceneView().getRootLayer();
        if (selectedLayer != null && rootLayer != selectedLayer) {
            moveLeft(selectedLayer);
        }
    }

    void moveLeft(Layer layer) {
        Layer parentLayer = layer.getParent();
        final Layer parentsParentLayer = parentLayer.getParent();
        if (parentsParentLayer != null ) {
            parentLayer.getChildren().remove(layer);
            final int parentIndex = parentsParentLayer.getChildIndex(parentLayer.getId());
            if (parentIndex < parentsParentLayer.getChildren().size() - 1) {
                parentsParentLayer.getChildren().add(parentIndex + 1, layer);
            } else {
                parentsParentLayer.getChildren().add(layer);
            }
        }
    }

}
