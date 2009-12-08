/*
 * $Id: ComputeRoiAreaAction.java,v 1.3 2009-12-08 16:08:33 lveci Exp $
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

import com.bc.ceres.core.ProgressMonitor;
import com.bc.ceres.swing.progress.DialogProgressMonitor;
import org.esa.beam.framework.datamodel.GeoCoding;
import org.esa.beam.framework.datamodel.GeoPos;
import org.esa.beam.framework.datamodel.Mask;
import org.esa.beam.framework.datamodel.PixelPos;
import org.esa.beam.framework.datamodel.Product;
import org.esa.beam.framework.datamodel.RasterDataNode;
import org.esa.beam.framework.ui.AbstractDialog;
import org.esa.beam.framework.ui.GridBagUtils;
import org.esa.beam.framework.ui.ModalDialog;
import org.esa.beam.framework.ui.command.CommandEvent;
import org.esa.beam.framework.ui.command.ExecCommand;
import org.esa.beam.framework.ui.product.ProductSceneView;
import org.esa.beam.util.Debug;
import org.esa.beam.util.math.MathUtils;
import org.esa.beam.util.math.RsMathUtils;
import org.esa.beam.visat.VisatApp;

import javax.swing.BoxLayout;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.SwingWorker;
import java.awt.Dialog;
import java.awt.GridBagConstraints;
import java.awt.Rectangle;
import java.awt.image.Raster;
import java.awt.image.RenderedImage;
import java.text.MessageFormat;
import java.util.concurrent.ExecutionException;

/**
 * Created by Marco Peters.
 *
 * @author Marco Peters
 * @version $Revision: 1.3 $ $Date: 2009-12-08 16:08:33 $
 */
public class ComputeRoiAreaAction extends ExecCommand {

    private static final String DIALOG_TITLE = "Compute ROI-Mask area"; /*I18N*/

    @Override
    public void actionPerformed(CommandEvent event) {
        computeROIArea();
    }

    @Override
    public void updateState(final CommandEvent event) {
        final ProductSceneView view = VisatApp.getApp().getSelectedProductSceneView();
        boolean enabled = false;

        if (view != null) {
            final RasterDataNode raster = view.getRaster();
            if (raster != null && raster.getGeoCoding() != null) {
                Product product = raster.getProduct();
                int numMasks = product.getMaskGroup().getNodeCount();
                enabled = numMasks > 0;
            }
        }
        setEnabled(enabled);
    }

    private void computeROIArea() {
        final String errMsgBase = "Failed to compute ROI-Mask area:\n";

        // Get selected VISAT view showing a product's band
        final ProductSceneView view = VisatApp.getApp().getSelectedProductSceneView();
        if (view == null) {
            VisatApp.getApp().showErrorDialog(DIALOG_TITLE, errMsgBase + "No view.");
            return;
        }

        // Get the current raster data node (band or tie-point grid)
        final RasterDataNode raster = view.getRaster();
        assert raster != null;
        
        String[] maskNames;
        if (raster.getRoiMaskGroup().getNodeCount() > 0) {
            maskNames = raster.getRoiMaskGroup().getNodeNames();
        } else {
            maskNames = raster.getProduct().getMaskGroup().getNodeNames();
        }
        String maskName;
        if (maskNames.length == 1) {
            maskName = maskNames[0];
        } else {
            JPanel panel = new JPanel();
            BoxLayout boxLayout = new BoxLayout(panel, BoxLayout.X_AXIS);
            panel.setLayout(boxLayout);
            panel.add(new JLabel("Select ROI-Mask: "));
            JComboBox maskCombo = new JComboBox(maskNames);
            panel.add(maskCombo);
            ModalDialog modalDialog = new ModalDialog(VisatApp.getApp().getApplicationWindow(), DIALOG_TITLE, panel, 
                                                      ModalDialog.ID_OK_CANCEL | ModalDialog.ID_HELP, getHelpId());
            if (modalDialog.show() == AbstractDialog.ID_OK) {
                maskName = (String) maskCombo.getSelectedItem();
            } else {
                return;
            }
        }
        Mask mask = raster.getProduct().getMaskGroup().get(maskName);
        
        RenderedImage maskImage = mask.getSourceImage();
        if (maskImage == null) {
            VisatApp.getApp().showErrorDialog(DIALOG_TITLE, errMsgBase + "No ROI-Mask image available.");
            return;
        }
            
        final SwingWorker<RoiAreaStatistics, Object> swingWorker = new RoiAreaSwingWorker(raster, maskImage, errMsgBase);
        swingWorker.execute();
    }


    private static class RoiAreaStatistics {

        private double earthRadius;
        private double roiArea;
        private double pixelAreaMin;
        private double pixelAreaMax;
        private int numPixels;

        private RoiAreaStatistics(double earthRadius) {
            this.earthRadius = earthRadius;
            roiArea = 0.0;
            pixelAreaMax = Double.NEGATIVE_INFINITY;
            pixelAreaMin  = Double.POSITIVE_INFINITY;
            numPixels = 0;
        }

        public double getEarthRadius() {
            return earthRadius;
        }

        public double getRoiArea() {
            return roiArea;
        }

        public void setRoiArea(double roiArea) {
            this.roiArea = roiArea;
        }

        public double getPixelAreaMin() {
            return pixelAreaMin;
        }

        public void setPixelAreaMin(double pixelAreaMin) {
            this.pixelAreaMin = pixelAreaMin;
        }

        public double getPixelAreaMax() {
            return pixelAreaMax;
        }

        public void setPixelAreaMax(double pixelAreaMax) {
            this.pixelAreaMax = pixelAreaMax;
        }

        public int getNumPixels() {
            return numPixels;
        }

        public void setNumPixels(int numPixels) {
            this.numPixels = numPixels;
        }
    }

    private class RoiAreaSwingWorker extends SwingWorker<RoiAreaStatistics, Object> {

        private final RasterDataNode raster;
        private final RenderedImage maskImage;
        private final String errMsgBase;

        private RoiAreaSwingWorker(RasterDataNode raster, RenderedImage maskImage, String errMsgBase) {
            this.raster = raster;
            this.maskImage = maskImage;
            this.errMsgBase = errMsgBase;
        }

        @Override
            protected RoiAreaStatistics doInBackground() throws Exception {
            ProgressMonitor pm = new DialogProgressMonitor(VisatApp.getApp().getMainFrame(), "Computing ROI-Mask area",
                                                           Dialog.ModalityType.APPLICATION_MODAL);
            return computeRoiAreaStatistics(pm);
        }

        private RoiAreaStatistics computeRoiAreaStatistics(ProgressMonitor pm) {
            
            final int minTileX = maskImage.getMinTileX();
            final int minTileY = maskImage.getMinTileY();

            final int numXTiles = maskImage.getNumXTiles();
            final int numYTiles = maskImage.getNumYTiles();
            
            final int w = raster.getSceneRasterWidth();
            final int h = raster.getSceneRasterHeight();
            final Rectangle imageRect = new Rectangle(0, 0, w, h);

            final PixelPos[] pixelPoints = new PixelPos[5];
            final GeoPos[] geoPoints = new GeoPos[5];
            for (int i = 0; i < geoPoints.length; i++) {
                pixelPoints[i] = new PixelPos();
                geoPoints[i] = new GeoPos();
            }

            RoiAreaStatistics areaStatistics = new RoiAreaStatistics(RsMathUtils.MEAN_EARTH_RADIUS / 1000.0);
            GeoCoding geoCoding = raster.getGeoCoding();

            pm.beginTask("Computing ROI-Mask area...", numXTiles * numYTiles);
            try {
                for (int tileX = minTileX; tileX < minTileX + numXTiles; ++tileX) {
                    for (int tileY = minTileY; tileY < minTileY + numYTiles; ++tileY) {
                        if (pm.isCanceled()) {
                            break;
                        }
                        final Rectangle tileRectangle = new Rectangle(
                                maskImage.getTileGridXOffset() + tileX * maskImage.getTileWidth(),
                                maskImage.getTileGridYOffset() + tileY * maskImage.getTileHeight(),
                                maskImage.getTileWidth(), maskImage.getTileHeight());
                        
                        final Rectangle r = imageRect.intersection(tileRectangle);
                        if (!r.isEmpty()) {
                            Raster roiTile = maskImage.getTile(tileX, tileY);
                            for (int y = r.y; y < r.y + r.height; y++) {
                                for (int x = r.x; x < r.x + r.width; x++) {
                                    if (roiTile.getSample(x, y, 0) != 0) {
                                        // 0: pixel center point
                                        pixelPoints[0].setLocation(x + 0.5f, y + 0.5f);
                                        // 1 --> 2 : parallel (geogr. hor. line) crossing pixel center point
                                        pixelPoints[1].setLocation(x + 0.0f, y + 0.5f);
                                        pixelPoints[2].setLocation(x + 1.0f, y + 0.5f);
                                        // 3 --> 4 : meridian (geogr. ver. line) crossing pixel center point
                                        pixelPoints[3].setLocation(x + 0.5f, y + 0.0f);
                                        pixelPoints[4].setLocation(x + 0.5f, y + 1.0f);

                                        for (int i = 0; i < geoPoints.length; i++) {
                                            geoCoding.getGeoPos(pixelPoints[i], geoPoints[i]);
                                        }
                                        float deltaLon = Math.abs(geoPoints[2].getLon() - geoPoints[1].getLon());
                                        float deltaLat = Math.abs(geoPoints[4].getLat() - geoPoints[3].getLat());
                                        double r2 = areaStatistics.getEarthRadius() * Math.cos(geoPoints[0].getLat() * MathUtils.DTOR);
                                        double a = r2 * deltaLon * MathUtils.DTOR;
                                        double b = areaStatistics.getEarthRadius() * deltaLat * MathUtils.DTOR;
                                        double pixelArea = a * b;
                                        areaStatistics.setPixelAreaMin(Math.min(areaStatistics.getPixelAreaMin(), pixelArea));
                                        areaStatistics.setPixelAreaMax(Math.max(areaStatistics.getPixelAreaMax(), pixelArea));
                                        areaStatistics.setRoiArea(areaStatistics.getRoiArea() + a * b);
                                        areaStatistics.setNumPixels(areaStatistics.getNumPixels() + 1);
                                    }
                                }
                            }
                        }
                        pm.worked(1);
                    }
                }
            } finally {
                pm.done();
            }
            return areaStatistics;
        }


        @Override
            public void done() {
            try {
                final RoiAreaStatistics areaStatistics = get();
                if (areaStatistics.getNumPixels() == 0) {
                    final String message = MessageFormat.format("{0}ROI-Mask is empty.", errMsgBase);
                    VisatApp.getApp().showErrorDialog(DIALOG_TITLE, message);
                } else {
                    showResults(areaStatistics);
                }
            } catch (ExecutionException e) {
                final String message = MessageFormat.format("An internal Error occured:\n{0}", e.getMessage());
                VisatApp.getApp().showErrorDialog(DIALOG_TITLE, message);
                Debug.trace(e);
            } catch (InterruptedException e) {
                final String message = MessageFormat.format("An internal Error occured:\n{0}", e.getMessage());
                VisatApp.getApp().showErrorDialog(DIALOG_TITLE, message);
                Debug.trace(e);
            }
        }

        private void showResults(RoiAreaStatistics areaStatistics) {
            final double roundFactor = 10000.0;
            final double roiAreaR = MathUtils.round(areaStatistics.getRoiArea(), roundFactor);
            final double meanPixelAreaR = MathUtils.round(areaStatistics.getRoiArea() / areaStatistics.getNumPixels(), roundFactor);
            final double pixelAreaMinR = MathUtils.round(areaStatistics.getPixelAreaMin(), roundFactor);
            final double pixelAreaMaxR = MathUtils.round(areaStatistics.getPixelAreaMax(), roundFactor);

            final JPanel content = GridBagUtils.createPanel();
            final GridBagConstraints gbc = new GridBagConstraints();
            gbc.fill = GridBagConstraints.HORIZONTAL;
            gbc.anchor = GridBagConstraints.WEST;
            gbc.insets.right = 4;
            gbc.gridy = 0;
            gbc.weightx = 0;

            gbc.insets.top = 2;
            addField(content, gbc, "Number of ROI pixels:", String.format("%15d",areaStatistics.getNumPixels()), "");
            addField(content, gbc, "ROI area:", String.format("%15.3f", roiAreaR), "km^2");
            addField(content, gbc, "Mean pixel area:", String.format("%15.3f", meanPixelAreaR), "km^2");
            addField(content, gbc, "Minimum pixel area:", String.format("%15.3f", pixelAreaMinR), "km^2");
            addField(content, gbc, "Maximum pixel area:", String.format("%15.3f", pixelAreaMaxR), "km^2");
            gbc.insets.top = 8;
            addField(content, gbc, "Mean earth radius:", String.format("%15.3f", areaStatistics.getEarthRadius()), "km");
            final ModalDialog dialog = new ModalDialog(VisatApp.getApp().getMainFrame(),
                                                       DIALOG_TITLE,
                                                       content,
                                                       ModalDialog.ID_OK | ModalDialog.ID_HELP,
                                                       getHelpId());
            dialog.show();

        }
        private void addField(final JPanel content, final GridBagConstraints gbc,
                                     final String text, final String value,
                                     final String unit) {
            content.add(new JLabel(text), gbc);
            gbc.weightx = 1;
            content.add(createTextField(value), gbc);
            gbc.weightx = 0;
            content.add(new JLabel(unit), gbc);
            gbc.gridy++;
        }

        private JTextField createTextField(final String value) {
            JTextField field = new JTextField(value);
            field.setEditable(false);
            field.setHorizontalAlignment(JTextField.RIGHT);
            return field;
        }


    }
}
