package org.esa.nest.dat.layers;

import com.bc.ceres.swing.figure.interactions.SelectionInteractor;
import com.bc.ceres.glayer.Layer;
import com.bc.ceres.glayer.LayerFilter;
import com.bc.ceres.glayer.support.LayerUtils;

import java.awt.event.MouseEvent;
import java.awt.*;
import java.util.List;

import org.esa.beam.framework.ui.product.ProductSceneView;
import org.esa.beam.visat.VisatApp;

public class SelectLayerInteractor extends SelectionInteractor {

    public SelectLayerInteractor() {

    }

    protected boolean isMouseOverSelection(MouseEvent event) {
        final boolean figureSelected = super.isMouseOverSelection(event);
        final boolean layerSelected = false;
        return figureSelected || layerSelected;
    }

    public SelectRectangleTool createSelectRectangleTool() {
        return new SelectLayerRectangleTool();        
    }

    public SelectPointTool createSelectPointTool() {
        return new SelectLayerPointTool();
    }

    private class SelectLayerPointTool extends SelectPointTool {
        @Override
        public void end(MouseEvent event) {
            super.end(event);

            final ProductSceneView view = VisatApp.getApp().getSelectedProductSceneView();
            final List<Layer> layers = findLayerSelections(view);

            for(Layer layer : layers) {
                LayerSelection laySel = (LayerSelection) layer;
                laySel.selectPoint(event.getX(), event.getY());
            }
            view.repaint();
        }
    }

    private class SelectLayerRectangleTool extends SelectRectangleTool {

        @Override
        public void drag(MouseEvent event) {
            super.drag(event);
            int width = event.getX() - referencePoint.x;
            int height = event.getY() - referencePoint.y;
            int x = referencePoint.x;
            int y = referencePoint.y;
            if (width < 0) {
                width *= -1;
                x -= width;
            }
            if (height < 0) {
                height *= -1;
                y -= height;
            }

            final ProductSceneView view = VisatApp.getApp().getSelectedProductSceneView();
            final List<Layer> layers = findLayerSelections(view);

            for(Layer layer : layers) {
                LayerSelection laySel = (LayerSelection) layer;
                laySel.selectRectangle(new Rectangle(x, y, width, height));
            }
        }

        @Override
        public void end(MouseEvent event) {
            super.end(event);
            final ProductSceneView view = VisatApp.getApp().getSelectedProductSceneView();
            view.repaint();
        }
    }

    private static List<Layer> findLayerSelections(final ProductSceneView view) {
        return LayerUtils.getChildLayers(view.getRootLayer(), LayerUtils.SearchMode.DEEP, new LayerFilter() {
            @Override
            public boolean accept(Layer layer) {
                return layer instanceof LayerSelection;
            }
        });
    }
}
