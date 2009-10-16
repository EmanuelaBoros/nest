package org.esa.beam.visat.toolviews.layermanager;

import com.bc.ceres.glayer.Layer;
import com.bc.ceres.glayer.support.LayerStyleListener;
import com.bc.ceres.glayer.support.LayerUtils;
import com.bc.ceres.swing.TreeCellExtender;
import com.jidesoft.swing.CheckBoxTree;
import com.jidesoft.swing.CheckBoxTreeSelectionModel;
import org.esa.beam.framework.ui.AppContext;
import org.esa.beam.framework.ui.UIUtils;
import org.esa.beam.framework.ui.GridBagUtils;
import org.esa.beam.framework.ui.product.ProductSceneView;
import org.esa.beam.framework.ui.tool.ToolButtonFactory;
import org.esa.beam.framework.help.HelpSys;
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
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.Rectangle;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.beans.PropertyChangeEvent;
import java.util.List;
import java.util.Hashtable;

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

    LayerManagerForm(AppContext appContext, String helpId) {
        super(appContext);
        this.view = appContext.getSelectedProductSceneView();
        initUI(helpId);
    }

    private void initUI(String helpId) {
        layerTreeModel = new LayerTreeModel(view.getRootLayer());
        layerTree = createCheckBoxTree(layerTreeModel);
        layerTree.setCellRenderer(new MyTreeCellRenderer());
        TreeCellExtender.equip(layerTree);

        Hashtable<Integer, JLabel> transparencySliderLabelTable = new Hashtable<Integer, JLabel>();
        transparencySliderLabelTable.put(0, createSliderLabel("0%"));
        transparencySliderLabelTable.put(127, createSliderLabel("50%"));
        transparencySliderLabelTable.put(255, createSliderLabel("100%"));
        transparencySlider = new JSlider(0, 255, 0);
        transparencySlider.setLabelTable(transparencySliderLabelTable);
        transparencySlider.setPaintLabels(true);
        transparencySlider.addChangeListener(new TransparencySliderListener());

        transparencyLabel = new JLabel("Transparency:");

        final JPanel sliderPanel = new JPanel(new BorderLayout(4, 4));
        sliderPanel.setBorder(new EmptyBorder(4, 4, 4, 4));
        sliderPanel.add(transparencyLabel, BorderLayout.WEST);
        sliderPanel.add(transparencySlider, BorderLayout.CENTER);

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

        AbstractButton helpButton = createToolButton("icons/Help24.gif");
        helpButton.setToolTipText("Help."); /*I18N*/
        helpButton.setName("helpButton");

        final JPanel actionBar = GridBagUtils.createPanel();
        final GridBagConstraints gbc = new GridBagConstraints();
        gbc.anchor = GridBagConstraints.CENTER;
        gbc.fill = GridBagConstraints.NONE;
        gbc.insets.top = 2;
        gbc.gridy = 0;
        actionBar.add(addButton, gbc);
        gbc.gridy++;
        actionBar.add(removeButton, gbc);
        gbc.gridy++;
        actionBar.add(openButton, gbc);
        gbc.gridy++;
        actionBar.add(upButton, gbc);
        gbc.gridy++;
        actionBar.add(downButton, gbc);
        gbc.gridy++;
        actionBar.add(leftButton, gbc);
        gbc.gridy++;
        actionBar.add(rightButton, gbc);
        gbc.gridy++;
        gbc.insets.bottom = 0;
        gbc.fill = GridBagConstraints.VERTICAL;
        gbc.weighty = 1.0;
        gbc.gridwidth = 2;
        actionBar.add(new JLabel(" "), gbc); // filler
        gbc.fill = GridBagConstraints.NONE;
        gbc.weighty = 0.0;
        gbc.gridy = 10;
        gbc.anchor = GridBagConstraints.EAST;
        actionBar.add(helpButton, gbc);


        JPanel layerPanel = new JPanel(new BorderLayout(4, 4));
        layerPanel.add(new JScrollPane(layerTree), BorderLayout.CENTER);
        layerPanel.add(sliderPanel, BorderLayout.SOUTH);

        control = new JPanel(new BorderLayout(4, 4));
        control.add(layerPanel, BorderLayout.CENTER);
        control.add(actionBar, BorderLayout.EAST);

        initLayerTreeVisibility(view.getRootLayer());
        updateFormControl();

        // todo - code duplication in all tool views with help support!!! (nf 200905)
        HelpSys.enableHelpOnButton(helpButton, helpId);
        HelpSys.enableHelpKey(control, helpId);

    }

    private static JLabel createSliderLabel(String text) {
        JLabel label = new JLabel(text);
        Font oldFont = label.getFont();
        Font newFont = oldFont.deriveFont(oldFont.getSize2D() * 0.85f);
        label.setFont(newFont);
        return label;
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

    public static boolean isLayerProtected(Layer layer) {
        return isLayerProtectedImpl(layer) || isChildLayerProtected(layer);
    }

    private Layer getSelectedLayer() {
        return view.getSelectedLayer();
    }

    private static boolean isLayerProtectedImpl(Layer layer) {
        return layer.getId().equals(ProductSceneView.BASE_IMAGE_LAYER_ID);
    }

    private static boolean isChildLayerProtected(Layer selectedLayer) {
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
            final double transparency = layer.getTransparency();
            final int n = (int) Math.round(255.0 * transparency);
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
                layer.setTransparency(transparencySlider.getValue() / 255.0);
                adjusting = false;
            }

        }
    }

    private class AddLayerActionListener implements ActionListener {
        Rectangle screenBounds;
        @Override
        public void actionPerformed(ActionEvent e) {
            LayerSourceAssistantPane pane = new LayerSourceAssistantPane(SwingUtilities.getWindowAncestor(control),
                                                                         "Add Layer",
                                                                         getAppContext());
            LayerSourceDescriptor[] layerSourceDescriptors = VisatActivator.getInstance().getLayerSources();
            pane.show(new SelectLayerSourceAssistantPage(layerSourceDescriptors), screenBounds);
            screenBounds = pane.getWindow().getBounds();
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

