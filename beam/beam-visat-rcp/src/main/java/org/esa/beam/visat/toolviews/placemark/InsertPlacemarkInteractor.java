package org.esa.beam.visat.toolviews.placemark;

import com.bc.ceres.swing.figure.FigureEditorInteractor;
import org.esa.beam.framework.datamodel.Pin;
import org.esa.beam.framework.datamodel.PixelPos;
import org.esa.beam.framework.datamodel.PlacemarkDescriptor;
import org.esa.beam.framework.datamodel.Product;
import org.esa.beam.framework.datamodel.PlacemarkGroup;
import org.esa.beam.framework.ui.product.ProductSceneView;

import java.awt.Cursor;
import java.awt.Image;
import java.awt.Toolkit;
import java.awt.event.MouseEvent;


public abstract class InsertPlacemarkInteractor extends FigureEditorInteractor {

    private final PlacemarkDescriptor placemarkDescriptor;
    private Cursor cursor;
    private boolean started;

    protected InsertPlacemarkInteractor(PlacemarkDescriptor placemarkDescriptor) {
        this.placemarkDescriptor = placemarkDescriptor;
        this.cursor = createCursor(placemarkDescriptor);
    }

    @Override
    public Cursor getCursor() {
        return cursor;
    }

    @Override
    public void mousePressed(MouseEvent event) {
        started = false;
        ProductSceneView sceneView = getProductSceneView(event);
        if (sceneView != null) {
            started = startInteraction(event);
        }
    }

    @Override
    public void mouseReleased(MouseEvent event) {
        if (started) {
            ProductSceneView sceneView = getProductSceneView(event);
            if (sceneView != null) {
                sceneView.selectVectorDataLayer(getPlacemarkGroup(sceneView.getProduct()).getVectorDataNode());
                if (isSingleButton1Click(event)) {
                    insertPlacemark(sceneView);
                }
                stopInteraction(event);
            }
        }
    }

    private void insertPlacemark(ProductSceneView view) {
        if (!view.isCurrentPixelPosValid()) {
            return;
        }
        Product product = view.getProduct();
        final String[] uniqueNameAndLabel = PlacemarkNameFactory.createUniqueNameAndLabel(placemarkDescriptor,
                                                                                          product);
        final String name = uniqueNameAndLabel[0];
        final String label = uniqueNameAndLabel[1];
        final PixelPos pixelPos = new PixelPos(view.getCurrentPixelX() + 0.5f,
                                               view.getCurrentPixelY() + 0.5f);
        final Pin newPlacemark = new Pin(name, label, "", pixelPos, null, placemarkDescriptor.createDefaultSymbol(),
                                         product.getGeoCoding());

        getPlacemarkGroup(product).add(newPlacemark);
    }

    private PlacemarkGroup getPlacemarkGroup(Product product) {
        return placemarkDescriptor.getPlacemarkGroup(product);
    }

    private static Cursor createCursor(PlacemarkDescriptor placemarkDescriptor) {
        final Image cursorImage = placemarkDescriptor.getCursorImage();
        if (cursorImage == null) {
            return Cursor.getPredefinedCursor(Cursor.CROSSHAIR_CURSOR);
        }
        return Toolkit.getDefaultToolkit().createCustomCursor(cursorImage,
                                                              placemarkDescriptor.getCursorHotSpot(),
                                                              placemarkDescriptor.getRoleName());
    }

    private ProductSceneView getProductSceneView(MouseEvent event) {
        if (event.getComponent() instanceof ProductSceneView) {
            return (ProductSceneView) event.getComponent();
        }
        if (event.getComponent().getParent() instanceof ProductSceneView) {
            return (ProductSceneView) event.getComponent().getParent();
        }
        return null;
    }
}