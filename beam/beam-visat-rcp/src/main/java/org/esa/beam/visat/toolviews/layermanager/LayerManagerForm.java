package org.esa.beam.visat.toolviews.layermanager;

import com.bc.ceres.glayer.Layer;
import com.bc.ceres.glayer.support.LayerStyleListener;
import com.bc.ceres.glayer.support.LayerUtils;
import com.jidesoft.swing.CheckBoxTree;
import com.jidesoft.swing.CheckBoxTreeSelectionModel;
import org.esa.beam.framework.ui.AppContext;
import org.esa.beam.framework.ui.UIUtils;
import org.esa.beam.framework.ui.product.ProductSceneView;
import org.esa.beam.framework.ui.tool.ToolButtonFactory;
import org.esa.beam.visat.VisatActivator;
import org.esa.beam.visat.toolviews.layermanager.layersrc.LayerSourceAssistantPane;
import org.esa.beam.visat.toolviews.layermanager.layersrc.SelectLayerSourceAssistantPage;

import javax.swing.AbstractButton;
import javax.swing.DropMode;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSlider;
import javax.swing.JTree;
import javax.swing.SwingUtilities;
import javax.swing.border.EmptyBorder;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.event.TreeSelectionEvent;
import javax.swing.event.TreeSelectionListener;
import javax.swing.tree.DefaultTreeCellRenderer;
import javax.swing.tree.TreePath;
import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.beans.PropertyChangeEvent;
import java.util.List;

class LayerManagerForm extends AbstractLayerForm {

    private final ProductSceneView view;
    private CheckBoxTree layerTree;
    private JSlider transparencySlider;
    private JPanel control;
    private boolean adjusting;
    private LayerTreeModel layerTreeModel;
    private JLabel transparencyLabel;
    private RemoveLayerAction removeLayerAction;
    private MoveLayerUpAction moveLayerUpAction;
    private MoveLayerDownAction moveLayerDownAction;
    private MoveLayerLeftAction moveLayerLeftAction;
    private MoveLayerRightAction moveLayerRightAction;
    private OpenLayerEditorAction openLayerEditorAction;

    LayerManagerForm(AppContext appContext) {
        super(appContext);
        this.view = appContext.getSelectedProductSceneView();
        initUI();
    }

    private void initUI() {
        layerTreeModel = new LayerTreeModel(view.getRootLayer());
        layerTree = createCheckBoxTree(layerTreeModel);
        layerTree.setCellRenderer(new MyTreeCellRenderer());

        transparencySlider = new JSlider(0, 100, 0);
        final JPanel sliderPanel = new JPanel(new BorderLayout(4, 4));
        sliderPanel.setBorder(new EmptyBorder(4, 4, 4, 4));
        transparencyLabel = new JLabel("Transparency:");
        sliderPanel.add(transparencyLabel, BorderLayout.WEST);
        sliderPanel.add(transparencySlider, BorderLayout.CENTER);
        transparencySlider.addChangeListener(new TransparencySliderListener());

        getRootLayer().addListener(new RootLayerListener());

        AbstractButton addButton = createToolButton("icons/Plus24.gif");
        addButton.addActionListener(new AddLayerActionListener());
        removeLayerAction = new RemoveLayerAction(getAppContext());
        AbstractButton removeButton = ToolButtonFactory.createButton(removeLayerAction, false);

        openLayerEditorAction = new OpenLayerEditorAction();
        AbstractButton openButton = ToolButtonFactory.createButton(openLayerEditorAction, false);

        moveLayerUpAction = new MoveLayerUpAction(getAppContext());
        AbstractButton upButton = ToolButtonFactory.createButton(moveLayerUpAction, false);

        moveLayerDownAction = new MoveLayerDownAction(getAppContext());
        AbstractButton downButton = ToolButtonFactory.createButton(moveLayerDownAction, false);

        moveLayerLeftAction = new MoveLayerLeftAction(getAppContext());
        AbstractButton leftButton = ToolButtonFactory.createButton(moveLayerLeftAction, false);

        moveLayerRightAction = new MoveLayerRightAction(getAppContext());
        AbstractButton rightButton = ToolButtonFactory.createButton(moveLayerRightAction, false);

        JPanel actionBar = new JPanel(new GridLayout(-1, 1, 2, 2));
        actionBar.add(addButton);
        actionBar.add(removeButton);
        actionBar.add(openButton);
        actionBar.add(upButton);
        actionBar.add(downButton);
        actionBar.add(leftButton);
        actionBar.add(rightButton);

        JPanel actionPanel = new JPanel(new BorderLayout());
        actionPanel.add(actionBar, BorderLayout.NORTH);

        control = new JPanel(new BorderLayout(4, 4));
        control.setBorder(new EmptyBorder(4, 4, 4, 4));
        control.add(new JScrollPane(layerTree), BorderLayout.CENTER);
        control.add(sliderPanel, BorderLayout.SOUTH);
        control.add(actionPanel, BorderLayout.EAST);

        initLayerTreeVisibility(view.getRootLayer());
        updateFormControl();
    }

    public Layer getRootLayer() {
        return view.getRootLayer();
    }

    @Override
    public JComponent getFormControl() {
        return control;
    }

    @Override
    public void updateFormControl() {
        Layer selectedLayer = getSelectedLayer();
        updateLayerStyleUI(selectedLayer);
        updateLayerTreeSelection(selectedLayer);
        boolean isLayerSelected = selectedLayer != null;
        removeLayerAction.setEnabled(isLayerSelected && !isLayerProtected(selectedLayer));
        openLayerEditorAction.setEnabled(isLayerSelected);
        moveLayerUpAction.setEnabled(isLayerSelected && moveLayerUpAction.canMove(selectedLayer));
        moveLayerDownAction.setEnabled(isLayerSelected && moveLayerDownAction.canMove(selectedLayer));
        moveLayerLeftAction.setEnabled(isLayerSelected && moveLayerLeftAction.canMove(selectedLayer));
        moveLayerRightAction.setEnabled(isLayerSelected && moveLayerRightAction.canMove(selectedLayer));
    }

    public boolean isLayerProtected(Layer layer) {
        return isLayerProtectedImpl(layer) || isChildLayerProtected(layer);
    }

    private Layer getSelectedLayer() {
        return view.getSelectedLayer();
    }

    private boolean isLayerProtectedImpl(Layer layer) {
        return layer.getId().equals(ProductSceneView.BASE_IMAGE_LAYER_ID);
    }

    private boolean isChildLayerProtected(Layer selectedLayer) {
        Layer[] children = selectedLayer.getChildren().toArray(new Layer[selectedLayer.getChildren().size()]);
        for (Layer childLayer : children) {
            if (isLayerProtectedImpl(childLayer) ||
                isChildLayerProtected(childLayer)) {
                return true;
            }
        }
        return false;
    }

    private static Layer getLayer(TreePath path) {
        if (path == null) {
            return null;
        }
        return (Layer) path.getLastPathComponent();
    }

    private void initLayerTreeVisibility(final Layer layer) {
        updateLayerTreeVisibility(layer);
        for (Layer childLayer : layer.getChildren()) {
            initLayerTreeVisibility(childLayer);
        }
    }

    private void updateLayerTreeVisibility(Layer layer) {
        CheckBoxTreeSelectionModel checkBoxTreeSelectionModel = layerTree.getCheckBoxTreeSelectionModel();
        Layer[] layerPath = LayerUtils.getLayerPath(layerTreeModel.getRootLayer(), layer);
        if (layerPath.length > 0) {
            if (layer.isVisible()) {
                checkBoxTreeSelectionModel.addSelectionPath(new TreePath(layerPath));
            } else {
                checkBoxTreeSelectionModel.removeSelectionPath(new TreePath(layerPath));
            }
            final List<Layer> children = layer.getChildren();
            if (!children.isEmpty()) {
                for (Layer child : children) {
                    updateLayerTreeVisibility(child);
                }
            }
        }
    }

    private void updateLayerTreeSelection(Layer selectedLayer) {
        if (selectedLayer != null) {
            Layer[] layerPath = LayerUtils.getLayerPath(layerTreeModel.getRootLayer(), selectedLayer);
            if (layerPath.length > 0) {
                layerTree.setSelectionPath(new TreePath(layerPath));
            } else {
                layerTree.clearSelection();
            }
        } else {
            layerTree.clearSelection();
        }
    }

    private void updateLayerStyleUI(Layer layer) {
        transparencyLabel.setEnabled(layer != null);
        transparencySlider.setEnabled(layer != null);
        if (layer != null) {
            final double transparency = 1 - layer.getStyle().getOpacity();
            final int n = (int) Math.round(100.0 * transparency);
            transparencySlider.setValue(n);
        }
    }

    private CheckBoxTree createCheckBoxTree(LayerTreeModel layerTreeModel) {

        final CheckBoxTree checkBoxTree = new CheckBoxTree(layerTreeModel);
        checkBoxTree.setRootVisible(false);
        checkBoxTree.setShowsRootHandles(true);
        checkBoxTree.setDigIn(false);

        checkBoxTree.setEditable(true);
        checkBoxTree.setDragEnabled(true);
        checkBoxTree.setDropMode(DropMode.ON_OR_INSERT);
        checkBoxTree.setTransferHandler(new LayerTreeTransferHandler(view, checkBoxTree));

        checkBoxTree.getSelectionModel().addTreeSelectionListener(new TreeSelectionListener() {
            @Override
            public void valueChanged(TreeSelectionEvent event) {
                Layer selectedLayer;
                final TreePath path = checkBoxTree.getSelectionPath();
                if (path != null) {
                    selectedLayer = getLayer(event.getPath());
                } else {
                    selectedLayer = null;
                }
                getAppContext().getSelectedProductSceneView().setSelectedLayer(selectedLayer);
            }
        });

        final CheckBoxTreeSelectionModel checkBoxSelectionModel = checkBoxTree.getCheckBoxTreeSelectionModel();
        checkBoxSelectionModel.addTreeSelectionListener(new TreeSelectionListener() {
            @Override
            public void valueChanged(TreeSelectionEvent event) {
                if (!adjusting) {
                    Layer layer = getLayer(event.getPath());
                    layer.setVisible(checkBoxSelectionModel.isPathSelected(event.getPath()));
                }
            }
        });

        final DefaultTreeCellRenderer renderer = (DefaultTreeCellRenderer) checkBoxTree.getActualCellRenderer();
        renderer.setLeafIcon(null);
        renderer.setClosedIcon(null);
        renderer.setOpenIcon(null);
        return checkBoxTree;
    }

    public static AbstractButton createToolButton(final String iconPath) {
        return ToolButtonFactory.createButton(UIUtils.loadImageIcon(iconPath), false);
    }

    private class RootLayerListener extends LayerStyleListener {

        @Override
        public void handleLayerStylePropertyChanged(Layer layer, PropertyChangeEvent event) {
            if (!adjusting) {
                updateFormControl();
            }
        }

        @Override
        public void handleLayerPropertyChanged(Layer layer, PropertyChangeEvent event) {
            if ("visible".equals(event.getPropertyName())) {
                updateLayerTreeVisibility(layer);
            }
        }

        @Override
        public void handleLayersAdded(Layer parentLayer, Layer[] childLayers) {
            for (Layer layer : childLayers) {
                updateLayerTreeVisibility(layer);
                updateLayerTreeSelection(layer);
            }
        }
    }

    private class TransparencySliderListener implements ChangeListener {

        @Override
        public void stateChanged(ChangeEvent e) {

            TreePath path = layerTree.getSelectionPath();
            if (path != null) {
                Layer layer = getLayer(path);
                adjusting = true;
                layer.getStyle().setOpacity(1.0 - transparencySlider.getValue() / 100.0f);
                adjusting = false;
            }

        }
    }

    private class AddLayerActionListener implements ActionListener {

        @Override
        public void actionPerformed(ActionEvent e) {
            LayerSourceAssistantPane pane = new LayerSourceAssistantPane(SwingUtilities.getWindowAncestor(control),
                                                                         "Add Layer",
                                                                         getAppContext());
            LayerSourceDescriptor[] layerSourceDescriptors = VisatActivator.getInstance().getLayerSources();
            pane.show(new SelectLayerSourceAssistantPage(layerSourceDescriptors));
        }
    }

    private static class MyTreeCellRenderer extends DefaultTreeCellRenderer {

        @Override
        public Component getTreeCellRendererComponent(JTree tree, Object value, boolean sel, boolean expanded,
                                                      boolean leaf,
                                                      int row, boolean hasFocus) {
            JLabel label = (JLabel) super.getTreeCellRendererComponent(tree,
                                                                       value, sel,
                                                                       expanded, leaf, row,
                                                                       hasFocus);
            Layer layer = (Layer) value;
            if (ProductSceneView.BASE_IMAGE_LAYER_ID.equals(layer.getId())) {
                label.setText(String.format("<html><b>%s</b></html>", layer.getName()));
            }
            return label;

        }
    }
}

