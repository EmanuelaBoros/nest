/*
 * $Id: ExportImageAction.java,v 1.8 2009-12-11 20:46:14 lveci Exp $
 *
 * Copyright (C) 2002 by Brockmann Consult (info@brockmann-consult.de)
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
package org.esa.beam.visat.actions;

import com.bc.ceres.binding.Property;
import com.bc.ceres.binding.PropertyContainer;
import com.bc.ceres.binding.PropertyDescriptor;
import com.bc.ceres.binding.accessors.DefaultPropertyAccessor;
import com.bc.ceres.binding.converters.IntegerConverter;
import com.bc.ceres.binding.swing.BindingContext;
import com.bc.ceres.glayer.support.ImageLayer;
import com.bc.ceres.grender.Viewport;
import com.bc.ceres.grender.support.BufferedImageRendering;
import com.bc.ceres.grender.support.DefaultViewport;
import org.esa.beam.framework.ui.PropertyPane;
import org.esa.beam.framework.ui.command.CommandEvent;
import org.esa.beam.framework.ui.product.ProductSceneView;
import org.esa.beam.util.io.BeamFileChooser;
import org.esa.beam.util.math.MathUtils;
import org.esa.beam.visat.VisatApp;

import javax.swing.BorderFactory;
import javax.swing.BoxLayout;
import javax.swing.ButtonGroup;
import javax.swing.JComponent;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import java.awt.Dimension;
import java.awt.Graphics2D;
import java.awt.GridLayout;
import java.awt.Rectangle;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.awt.image.RenderedImage;

/**
 * Action for exporting scene views as images.
 *
 * @author Marco Peters
 * @author Ralf Quast
 * @version $Revision: 1.8 $ $Date: 2009-12-11 20:46:14 $
 */
public class ExportImageAction extends AbstractExportImageAction {

    private JRadioButton buttonFullScene;
    private SizeComponent sizeComponent;

    @Override
    public void actionPerformed(CommandEvent event) {
        exportImage(getVisatApp(), getSceneImageFileFilters(), event.getSelectableCommand());
    }

    @Override
    public void updateState(final CommandEvent event) {
        boolean enabled = getVisatApp().getSelectedProductSceneView() != null;
        event.getSelectableCommand().setEnabled(enabled);

    }

    @Override
    protected void configureFileChooser(final BeamFileChooser fileChooser, final ProductSceneView view,
                                        String imageBaseName) {
        fileChooser.setDialogTitle(getVisatApp().getAppName() + " - " + "Export Image"); /*I18N*/
        if (view.isRGB()) {
            fileChooser.setCurrentFilename(imageBaseName + "_RGB");
        } else {
            fileChooser.setCurrentFilename(imageBaseName + "_" + view.getRaster().getName());
        }
        final JPanel regionPanel = new JPanel(new GridLayout(2, 1));
        regionPanel.setBorder(BorderFactory.createTitledBorder("Image Region")); /*I18N*/
        buttonFullScene = new JRadioButton("Full scene", false);
        final JRadioButton buttonVisibleRegion = new JRadioButton("Visible region", true); /*I18N*/
        ButtonGroup buttonGroup = new ButtonGroup();
        buttonGroup.add(buttonVisibleRegion);
        buttonGroup.add(buttonFullScene);
        regionPanel.add(buttonVisibleRegion);
        regionPanel.add(buttonFullScene);
        sizeComponent = new SizeComponent(view);
        JComponent sizePanel = sizeComponent.createComponent();
        sizePanel.setBorder(BorderFactory.createTitledBorder("Image Dimension")); /*I18N*/
        final JPanel accessory = new JPanel();
        accessory.setLayout(new BoxLayout(accessory, BoxLayout.Y_AXIS));
        accessory.add(regionPanel);
        accessory.add(sizePanel);
        fileChooser.setAccessory(accessory);

        buttonFullScene.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                sizeComponent.updateDimensions();
            }
        });
        buttonVisibleRegion.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                sizeComponent.updateDimensions();
            }
        });

    }

    @Override
    protected RenderedImage createImage(String imageFormat, ProductSceneView view) {
        final boolean useAlpha = !"BMP".equals(imageFormat) && !"JPEG".equals(imageFormat);
        final boolean entireImage = isEntireImageSelected();

        return createImage(view, entireImage, sizeComponent.getDimension(), useAlpha);
    }

    static RenderedImage createImage(ProductSceneView view, boolean fullScene, Dimension dimension,
                                     boolean alphaChannel) {
        final int imageType = alphaChannel ? BufferedImage.TYPE_4BYTE_ABGR : BufferedImage.TYPE_3BYTE_BGR;
        final BufferedImage bufferedImage = new BufferedImage(dimension.width, dimension.height, imageType);

        final Viewport vp1 = view.getLayerCanvas().getViewport();
        final Viewport vp2 = new DefaultViewport(new Rectangle(dimension), vp1.isModelYAxisDown());
        if (fullScene) {
            vp2.zoom(view.getBaseImageLayer().getModelBounds());
        } else {
            setTransform(vp1, vp2);
        }

        final BufferedImageRendering imageRendering = new BufferedImageRendering(bufferedImage, vp2);
        if (!alphaChannel) {
            final Graphics2D graphics = imageRendering.getGraphics();
            graphics.setColor(view.getLayerCanvas().getBackground());
            graphics.fillRect(0, 0, dimension.width, dimension.height);
        }
        view.getRootLayer().render(imageRendering);

        return bufferedImage;
    }

    private static void setTransform(Viewport vp1, Viewport vp2) {
        vp2.setTransform(vp1);

        final Rectangle rectangle1 = vp1.getViewBounds();
        final Rectangle rectangle2 = vp2.getViewBounds();

        final double w1 = rectangle1.getWidth();
        final double w2 = rectangle2.getWidth();
        final double h1 = rectangle1.getHeight();
        final double h2 = rectangle2.getHeight();
        final double x1 = rectangle1.getX();
        final double y1 = rectangle1.getY();
        final double cx = (x1 + w1) / 2.0;
        final double cy = (y1 + h1) / 2.0;

        final double magnification;
        if (w1 > h1) {
            magnification = w2 / w1;
        } else {
            magnification = h2 / h1;
        }

        final Point2D modelCenter = vp1.getViewToModelTransform().transform(new Point2D.Double(cx, cy), null);
        vp2.setZoomFactor(vp1.getZoomFactor() * magnification, modelCenter.getX(), modelCenter.getY());
    }

    @Override
    protected boolean isEntireImageSelected() {
        return buttonFullScene.isSelected();
    }

    private class SizeComponent {

        private static final String PROPERTY_NAME_HEIGHT = "height";
        private static final String PROPERTY_NAME_WIDTH = "width";

        private final PropertyContainer propertyContainer;
        private final ProductSceneView view;

        public SizeComponent(ProductSceneView view) {
            this.view = view;
            propertyContainer = new PropertyContainer();
            initValueContainer();
            updateDimensions();
        }

        public void updateDimensions() {
            final Rectangle2D bounds;
            if (isEntireImageSelected()) {
                final ImageLayer imageLayer = view.getBaseImageLayer();
                final Rectangle2D modelBounds = imageLayer.getModelBounds();
                bounds = imageLayer.getModelToImageTransform().createTransformedShape(modelBounds).getBounds2D();
            } else {
                bounds = view.getLayerCanvas().getViewport().getViewBounds();
            }

            int w = toInteger(bounds.getWidth());
            int h = toInteger(bounds.getHeight());

            final long freeMemory = getFreeMemory();
            final long expectedMemory = getExpectedMemory(w, h);
            if (freeMemory < expectedMemory) {
                final int answer = showQuestionDialog();
                if (answer != JOptionPane.YES_OPTION) {
                    final double scale = Math.sqrt((double) freeMemory / (double) expectedMemory);
                    final double scaledW = w * scale;
                    final double scaledH = h * scale;

                    w = toInteger(scaledW);
                    h = toInteger(scaledH);
                }

            }

            setWidth(w);
            setHeight(h);
        }

        private int toInteger(double value) {
            return MathUtils.floorInt(value);
        }

        public JComponent createComponent() {
            BindingContext bindingContext = new BindingContext(propertyContainer);
            PropertyPane propertyPane = new PropertyPane(bindingContext);
            return propertyPane.createPanel();
        }

        public Dimension getDimension() {
            return new Dimension(getWidth(), getHeight());
        }

        private void initValueContainer() {
            final PropertyDescriptor widthDescriptor = new PropertyDescriptor(PROPERTY_NAME_WIDTH, Integer.class);
            widthDescriptor.setConverter(new IntegerConverter());
            propertyContainer.addProperty(new Property(widthDescriptor, new DefaultPropertyAccessor()));

            final PropertyDescriptor heightDescriptor = new PropertyDescriptor(PROPERTY_NAME_HEIGHT, Integer.class);
            heightDescriptor.setConverter(new IntegerConverter());
            propertyContainer.addProperty(new Property(heightDescriptor, new DefaultPropertyAccessor()));
        }

        private int showQuestionDialog() {
            return VisatApp.getApp().showQuestionDialog(
                    "There may not be enough memory to export the image because\n" +
                    "the image dimension is too large.\n\n" +
                    "Do you really want to keep the image dimension?", null);
        }

        private long getFreeMemory() {
            final long usedMemory = Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory();
            return Runtime.getRuntime().maxMemory() - usedMemory;
        }

        private long getExpectedMemory(int width, int height) {
            return width * height * 6L;
        }

        private int getWidth() {
            return (Integer) propertyContainer.getValue(PROPERTY_NAME_WIDTH);
        }

        private void setWidth(Object value) {
            propertyContainer.setValue(PROPERTY_NAME_WIDTH, value);
        }

        private int getHeight() {
            return (Integer) propertyContainer.getValue(PROPERTY_NAME_HEIGHT);
        }

        private void setHeight(Object value) {
            propertyContainer.setValue(PROPERTY_NAME_HEIGHT, value);
        }
    }
}
