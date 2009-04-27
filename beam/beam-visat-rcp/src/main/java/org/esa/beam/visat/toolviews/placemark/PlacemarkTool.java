package org.esa.beam.visat.toolviews.placemark;

import org.esa.beam.framework.datamodel.*;
import org.esa.beam.framework.draw.Drawable;
import org.esa.beam.framework.ui.UIUtils;
import org.esa.beam.framework.ui.command.Command;
import org.esa.beam.framework.ui.command.ExecCommand;
import org.esa.beam.framework.ui.product.ProductSceneView;
import org.esa.beam.framework.ui.tool.AbstractTool;
import org.esa.beam.framework.ui.tool.DrawingEditor;
import org.esa.beam.framework.ui.tool.ToolInputEvent;
import org.esa.beam.visat.VisatApp;

import javax.swing.*;
import java.awt.*;
import java.awt.geom.AffineTransform;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.util.Collection;

/**
 * Created by Marco Peters.
 *
 * @author Marco Peters
 * @version $Revision: 1.1 $ $Date: 2009-04-27 13:08:26 $
 */
public abstract class PlacemarkTool extends AbstractTool {

    private final PlacemarkDescriptor placemarkDescriptor;
    private Pin draggedPlacemark;
    private Cursor cursor;

    protected PlacemarkTool(PlacemarkDescriptor placemarkDescriptor) {
        this.placemarkDescriptor = placemarkDescriptor;
        cursor = createCursor(placemarkDescriptor);
    }

    /**
     * Gets a thing that can be drawn while the tool is working.
     *
     * @return always <code>null</code>
     */
    @Override
    public Drawable getDrawable() {
        return null;
    }

    @Override
    public void mouseClicked(ToolInputEvent e) {
        activateOverlay();
        if (isSingleLeftClick(e)) {
            selectOrInsertPlacemark(e);
        } else if (isDoubleLeftClick(e)) {
            selectAndEditPlacemark(e);
        } else if (e.getMouseEvent().isPopupTrigger()) {
            showPopupMenu(e);
        }
    }

    @Override
    public void mousePressed(ToolInputEvent e) {
        Pin draggedPlacemark = getPlacemarkForPixelPos(e.getPixelX(), e.getPixelY());
        setDraggedPlacemark(draggedPlacemark);
    }


    @Override
    public void mouseDragged(ToolInputEvent e) {
        if (getDraggedPlacemark() != null && e.isPixelPosValid()) {
            PixelPos pixelPos = new PixelPos((float) e.getPixelPos().getX() + 0.5f,
                                             (float) e.getPixelPos().getY() + 0.5f);
            if (pixelPos.isValid()) {
                getDraggedPlacemark().setPixelPos(pixelPos);
            }
        }
    }

    private void selectOrInsertPlacemark(ToolInputEvent e) {
        if (!e.isPixelPosValid()) {
            return;
        }
        ProductSceneView view = getProductSceneView();
        if (view == null) {
            return;
        }
        Product product = view.getProduct();
        Pin clickedPlacemark = getPlacemarkForPixelPos(e.getPixelX(), e.getPixelY());
        if (clickedPlacemark != null) {
            setPlacemarkSelected(getPlacemarkGroup(product), clickedPlacemark, false);
        } else {
            final String[] uniqueNameAndLabel = PlacemarkNameFactory.createUniqueNameAndLabel(placemarkDescriptor,
                                                                                              product);
            final String name = uniqueNameAndLabel[0];
            final String label = uniqueNameAndLabel[1];
            final Pin newPlacemark = new Pin(name, label, "",
                                             new PixelPos(0.5f + e.getPixelX(), 0.5f + e.getPixelY()),
                                             null,
                                             placemarkDescriptor.createDefaultSymbol());
            getPlacemarkGroup(product).add(newPlacemark);
        }
    }

    private void selectAndEditPlacemark(ToolInputEvent e) {
        if (!e.isPixelPosValid()) {
            return;
        }
        ProductSceneView view = getProductSceneView();
        if (view == null) {
            return;
        }
        final Product product = view.getProduct();
        int pixelX = e.getPixelX();
        int pixelY = e.getPixelY();
        Pin clickedPlacemark = getPlacemarkForPixelPos(pixelX, pixelY);
        if (clickedPlacemark != null) {
            setPlacemarkSelected(getPlacemarkGroup(product), clickedPlacemark, true);
            boolean ok = PlacemarkDialog.showEditPlacemarkDialog(VisatApp.getApp().getMainFrame(), product,
                                                     clickedPlacemark, placemarkDescriptor);
            if (ok) {
                updateState();
            }
        }
    }

    private ExecCommand getShowOverlayCommand() {
        Command command = VisatApp.getApp().getCommandManager().getCommand(getShowOverlayCommandId());
        if (command != null) {
            if (command instanceof ExecCommand) {
                return (ExecCommand) command;
            }
        }
        return null;
    }

    protected String getShowOverlayCommandId() {
        return placemarkDescriptor.getShowLayerCommandId();
    }

    protected ProductNodeGroup<Pin> getPlacemarkGroup(Product product) {
        return placemarkDescriptor.getPlacemarkGroup(product);
    }

    @Override
    public Cursor getCursor() {
        return cursor;
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

    protected void showPopupMenu(ToolInputEvent e) {
        JPopupMenu popup = new JPopupMenu();
        VisatApp.getApp().getCommandUIFactory().addContextDependentMenuItems("pin", popup);
        UIUtils.showPopup(popup, e.getMouseEvent());
    }

    protected ProductSceneView getProductSceneView() {
        DrawingEditor drawingEditor = getDrawingEditor();
        if (drawingEditor instanceof ProductSceneView) {
            return (ProductSceneView) drawingEditor;
        }
        return null;
    }

    public static void setPlacemarkSelected(ProductNodeGroup<Pin> placemarkGroup, Pin clickedPlacemark,
                                            boolean forceSelection) {
        boolean select = true;
        Collection<Pin> selectedPlacemark = placemarkGroup.getSelectedNodes();
        for (Pin placemark : selectedPlacemark) {
            if (placemark.isSelected()) {
                if (placemark == clickedPlacemark) {
                    select = false;
                }
                placemark.setSelected(false);
            }
        }
        if (forceSelection || select) {
            clickedPlacemark.setSelected(true);
        }
        updateState();
    }

    private static void updateState() {
        VisatApp.getApp().updateState();
    }

    private void activateOverlay() {
        ExecCommand overlayCommand = getShowOverlayCommand();
        if (overlayCommand != null) {
            overlayCommand.setSelected(true);
            overlayCommand.execute();
        }
    }

    private Pin getPlacemarkForPixelPos(final int selPixelX, final int selPixelY) {
        ProductSceneView view = getProductSceneView();

        final AffineTransform i2v = view.getBaseImageToViewTransform();

        final Rectangle2D selViewRect = i2v.createTransformedShape(new Rectangle(selPixelX, selPixelY, 1, 1)).getBounds2D();

        // Compare with pin insertion points (which are in product raster pixel coordinates)
        //
        ProductNodeGroup<Pin> placemarkGroup = getPlacemarkGroup(getProductSceneView().getProduct());
        Pin[] placemarks = placemarkGroup.toArray(new Pin[placemarkGroup.getNodeCount()]);
        for (final Pin placemark : placemarks) {
            // Convert pin pixel to view coordinates
            PixelPos placemarkPixelPos = placemark.getPixelPos();
            if (placemarkPixelPos == null) {
                return null;
            }
            final Rectangle2D placemarkViewRect = i2v.createTransformedShape(new Rectangle2D.Double(Math.floor(placemarkPixelPos.getX()),
                                                                                                    Math.floor(placemarkPixelPos.getY()),
                                                                                                    1, 1)).getBounds2D();
            // Use a rectangular region around the insertion point for comparision
            if (selViewRect.intersects(placemarkViewRect)) {
                return placemark;
            }
        }

        // Now compare against pin symbols (which are in view coordinates).
        // Depending on the symbol used, this may be more or less expensive.
        //
        for (final Pin placemark : placemarks) {
            PlacemarkSymbol symbol = placemark.getSymbol();

            // Convert pin pixel to view coordinates
            PixelPos placemarkPixelPos = placemark.getPixelPos();
            final Point2D placemarkViewPos = i2v.transform(new Point2D.Double(placemarkPixelPos.getX(),
                                                                              placemarkPixelPos.getY()),
                                                           null);

            PixelPos refPoint = symbol.getRefPoint();
            if (refPoint == null) {
                refPoint = new PixelPos(0, 0);
            }
            final Rectangle2D.Double relViewRect = new Rectangle2D.Double(selViewRect.getX() - placemarkViewPos.getX() + refPoint.getX(),
                                                                          selViewRect.getY() - placemarkViewPos.getY() + refPoint.getY(),
                                                                          selViewRect.getWidth(),
                                                                          selViewRect.getHeight());

            if (symbol.getShape().intersects(relViewRect)) {
                return placemark;
            }
        }
        return null;
    }

    protected Pin getDraggedPlacemark() {
        return draggedPlacemark;
    }

    protected void setDraggedPlacemark(Pin draggedPlacemark) {
        this.draggedPlacemark = draggedPlacemark;
    }
}
