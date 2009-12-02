package org.esa.beam.visat.actions;

import com.bc.ceres.core.ProgressMonitor;
import com.bc.ceres.swing.progress.ProgressMonitorSwingWorker;
import org.esa.beam.framework.datamodel.Band;
import org.esa.beam.framework.datamodel.GeoCoding;
import org.esa.beam.framework.datamodel.GeoPos;
import org.esa.beam.framework.datamodel.Mask;
import org.esa.beam.framework.datamodel.PixelPos;
import org.esa.beam.framework.datamodel.Product;
import org.esa.beam.framework.datamodel.RasterDataNode;
import org.esa.beam.framework.datamodel.TiePointGrid;
import org.esa.beam.framework.ui.AbstractDialog;
import org.esa.beam.framework.ui.ModalDialog;
import org.esa.beam.framework.ui.SelectExportMethodDialog;
import org.esa.beam.framework.ui.UIUtils;
import org.esa.beam.framework.ui.command.CommandEvent;
import org.esa.beam.framework.ui.command.ExecCommand;
import org.esa.beam.framework.ui.product.ProductSceneView;
import org.esa.beam.util.SystemUtils;
import org.esa.beam.util.io.FileUtils;
import org.esa.beam.visat.VisatApp;

import java.awt.Rectangle;
import java.awt.image.Raster;
import java.awt.image.RenderedImage;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;

import javax.swing.BoxLayout;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;

public class ExportROIPixelsAction extends ExecCommand {

    private final static String DLG_TITLE = "Export ROI-Mask Pixels";
    private static final String ERR_MSG_BASE = "ROI-Mask pixels cannot be exported:\n";

    /**
     * Invoked when a command action is performed.
     *
     * @param event the command event
     */
    @Override
    public void actionPerformed(CommandEvent event) {
        exportROIPixels();
    }

    /**
     * Called when a command should update its state.
     * <p/>
     * <p> This method can contain some code which analyzes the underlying element and makes a decision whether
     * this item or group should be made visible/invisible or enabled/disabled etc.
     *
     * @param event the command event
     */
    @Override
    public void updateState(final CommandEvent event) {
        final boolean enabled = hasSelectedRasterMasks();
        setEnabled(enabled);
    }
    /////////////////////////////////////////////////////////////////////////
    // Private implementations for the "Export ROI Pixels" command
    /////////////////////////////////////////////////////////////////////////

    /**
     * Performs the actual "Export ROI Pixels" command.
     */
    private void exportROIPixels() {

        if (!hasSelectedRasterMasks()) {
            VisatApp.getApp().showErrorDialog(DLG_TITLE,
                                              ERR_MSG_BASE + "There are no masks available in the currently selected product.");  /*I18N*/
            return;
        }

        // Get current VISAT view showing a product's band
        final ProductSceneView view = VisatApp.getApp().getSelectedProductSceneView();
        if (view == null) {
            return;
        }
        // Get the displayed raster data node (band or tie-point grid)
        final RasterDataNode raster = view.getRaster();
        
        String[] maskNames = raster.getProduct().getMaskGroup().getNodeNames();
        JPanel panel = new JPanel();
        BoxLayout boxLayout = new BoxLayout(panel, BoxLayout.X_AXIS);
        panel.setLayout(boxLayout);
        panel.add(new JLabel("Select ROI-Mask: "));
        JComboBox maskCombo = new JComboBox(maskNames);
        panel.add(maskCombo);
        ModalDialog modalDialog = new ModalDialog(VisatApp.getApp().getApplicationWindow(), DLG_TITLE, panel, 
                        ModalDialog.ID_OK_CANCEL | ModalDialog.ID_HELP, getHelpId());
        if (modalDialog.show() != AbstractDialog.ID_OK) {
            return; //Cancel
        }
        String maskName = (String) maskCombo.getSelectedItem();
        Mask mask = raster.getProduct().getMaskGroup().get(maskName);
        if (mask == null) {
            VisatApp.getApp().showErrorDialog(DLG_TITLE, ERR_MSG_BASE + "No ROI-Mask selected.");
            return;
        }
        final RenderedImage roiMaskImage = mask.getSourceImage();
        if (roiMaskImage == null) {
            VisatApp.getApp().showErrorDialog(DLG_TITLE, ERR_MSG_BASE + "No ROI-Mask image available.");
            return;
        }
        // Compute total number of ROI pixels
        final long numROIPixels = getNumROIPixels(raster, roiMaskImage);
        
        final String questionText = "How do you want to export the pixel values?\n"; /*I18N*/
        String numPixelsText;
        if (numROIPixels == 1) {
            numPixelsText = "One ROI-Mask pixel will be exported.\n"; /*I18N*/
        } else {
            numPixelsText = numROIPixels + " ROI-Mask pixels will be exported.\n"; /*I18N*/
        }
        // Get export method from user
        final int method = SelectExportMethodDialog.run(VisatApp.getApp().getMainFrame(), getWindowTitle(),
                                                        questionText + numPixelsText, getHelpId());
        final PrintWriter out;
        final StringBuffer clipboardText;
        final int initialBufferSize = 256000;
        if (method == SelectExportMethodDialog.EXPORT_TO_CLIPBOARD) {
            // Write into string buffer
            final StringWriter stringWriter = new StringWriter(initialBufferSize);
            out = new PrintWriter(stringWriter);
            clipboardText = stringWriter.getBuffer();
        } else if (method == SelectExportMethodDialog.EXPORT_TO_FILE) {
            // Write into file, get file from user
            final File file = promptForFile(VisatApp.getApp(), createDefaultFileName(raster, maskName));
            if (file == null) {
                return; // Cancel
            }
            final FileWriter fileWriter;
            try {
                fileWriter = new FileWriter(file);
            } catch (IOException e) {
                VisatApp.getApp().showErrorDialog(DLG_TITLE,
                                                  ERR_MSG_BASE + "Failed to create file '" + file + "':\n" + e.getMessage()); /*I18N*/
                return; // Error
            }
            out = new PrintWriter(new BufferedWriter(fileWriter, initialBufferSize));
            clipboardText = null;
        } else {
            return; // Cancel
        }

        final ProgressMonitorSwingWorker<Exception, Object> swingWorker = new ProgressMonitorSwingWorker<Exception, Object>(
                VisatApp.getApp().getMainFrame(), DLG_TITLE) {

            @Override
            protected Exception doInBackground(ProgressMonitor pm) throws Exception {
                boolean success;
                Exception returnValue = null;
                try {
                    success = exportROIPixels(out, raster.getProduct(), roiMaskImage, pm);
                    if (success && clipboardText != null) {
                        SystemUtils.copyToClipboard(clipboardText.toString());
                        clipboardText.setLength(0);
                    }
                } catch (Exception e) {
                    returnValue = e;
                } finally {
                    out.close();
                }
                return returnValue;
            }

            /**
             * Called on the event dispatching thread (not on the worker thread) after the <code>construct</code> method
             * has returned.
             */
            @Override
            public void done() {
                // clear status bar
                VisatApp.getApp().clearStatusBarMessage();
                // show default-cursor
                UIUtils.setRootFrameDefaultCursor(VisatApp.getApp().getMainFrame());
                // On error, show error message
                Exception exception;
                try {
                    exception = get();
                } catch (Exception e) {
                    exception = e;
                }
                if (exception != null) {
                    VisatApp.getApp().showErrorDialog(DLG_TITLE,
                                                      ERR_MSG_BASE + exception.getMessage());
                }
            }

        };

        // show wait-cursor
        UIUtils.setRootFrameWaitCursor(VisatApp.getApp().getMainFrame());
        // show message in status bar
        VisatApp.getApp().setStatusBarMessage("Exporting ROI pixels..."); /*I18N*/

        // Start separate worker thread.
        swingWorker.execute();
    }

    private static String createDefaultFileName(final RasterDataNode raster, String maskName) {
        String productName = FileUtils.getFilenameWithoutExtension(raster.getProduct().getName());
        String rasterName = raster.getName();
        return productName + "_" + rasterName + "_" + maskName + "_ROI.txt";
    }

    private static String getWindowTitle() {
        return VisatApp.getApp().getAppName() + " - " + DLG_TITLE;
    }

    /**
     * Opens a modal file chooser dialog that prompts the user to select the output file name.
     *
     * @param visatApp the VISAT application
     * @return the selected file, <code>null</code> means "Cancel"
     */
    private static File promptForFile(final VisatApp visatApp, String defaultFileName) {
        return visatApp.showFileSaveDialog(DLG_TITLE,
                                           false,
                                           null,
                                           ".txt",
                                           defaultFileName,
                                           "exportROIPixels.lastDir");
    }

    /**
     * Writes all pixel values of the given product within the given ROI to the specified out.
     *
     * @param out      the data output writer
     * @param product  the product providing the pixel values
     * @param roiImage the mask image for the ROI
     * @return <code>true</code> for success, <code>false</code> if export has been terminated (by user)
     */
    private static boolean exportROIPixels(final PrintWriter out,
                                           final Product product,
                                           final RenderedImage roiImage, ProgressMonitor pm) throws IOException {

        final Band[] bands = product.getBands();
        final TiePointGrid[] tiePointGrids = product.getTiePointGrids();
        final GeoCoding geoCoding = product.getGeoCoding();
        
        final int minTileX = roiImage.getMinTileX();
        final int minTileY = roiImage.getMinTileY();
        
        final int numXTiles = roiImage.getNumXTiles();
        final int numYTiles = roiImage.getNumYTiles();
        
        final int w = product.getSceneRasterWidth();
        final int h = product.getSceneRasterHeight();
        final Rectangle imageRect = new Rectangle(0, 0, w, h);
        
        pm.beginTask("Writing pixel data...", numXTiles * numYTiles + 1);
        try {
            writeHeaderLine(out, geoCoding, bands, tiePointGrids);
            pm.worked(1);

            for (int tileX = minTileX; tileX < minTileX + numXTiles; ++tileX) {
                for (int tileY = minTileY; tileY < minTileY + numYTiles; ++tileY) {
                    if (pm.isCanceled()) {
                        return false;
                    }
                    final Rectangle tileRectangle = new Rectangle(roiImage.getTileGridXOffset() + tileX * roiImage.getTileWidth(),
                                                                  roiImage.getTileGridYOffset() + tileY * roiImage.getTileHeight(),
                                                                  roiImage.getTileWidth(), roiImage.getTileHeight());

                    final Rectangle r = imageRect.intersection(tileRectangle);
                    if (!r.isEmpty()) {
                        Raster roiTile = roiImage.getTile(tileX, tileY);
                        for (int y = r.y; y < r.y + r.height; y++) {
                            for (int x = r.x; x < r.x + r.width; x++) {
                                if (roiTile.getSample(x, y, 0) != 0) {
                                    writeDataLine(out, geoCoding, bands, tiePointGrids, x, y);
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

        return true;
    }

    /**
     * Writes the header line of the dataset to be exported.
     *
     * @param out           the data output writer
     * @param geoCoding     the product's geo-coding
     * @param bands         the array of bands to be considered
     * @param tiePointGrids
     */
    private static void writeHeaderLine(final PrintWriter out,
                                        final GeoCoding geoCoding,
                                        final Band[] bands, TiePointGrid[] tiePointGrids) {
        out.print("Pixel-X");
        out.print("\t");
        out.print("Pixel-Y");
        if (geoCoding != null) {
            out.print("\t");
            out.print("Longitude");
            out.print("\t");
            out.print("Latitude");
        }
        for (final Band band : bands) {
            out.print("\t");
            out.print(band.getName());
        }
        for (final TiePointGrid tiePointGrid : tiePointGrids) {
            out.print("\t");
            out.print(tiePointGrid.getName());
        }
        out.print("\n");
    }

    /**
     * Writes a data line of the dataset to be exported for the given pixel position.
     *
     * @param out           the data output writer
     * @param geoCoding     the product's geo-coding
     * @param bands         the array of bands that provide pixel values
     * @param tiePointGrids
     * @param x             the current pixel's X coordinate
     * @param y             the current pixel's Y coordinate
     */
    private static void writeDataLine(final PrintWriter out,
                                      final GeoCoding geoCoding,
                                      final Band[] bands,
                                      TiePointGrid[] tiePointGrids, int x,
                                      int y) throws IOException {
        final PixelPos pixelPos = new PixelPos(x + 0.5f, y + 0.5f);
        final int[] intPixel = new int[1];
        final float[] floatPixel = new float[1];

        out.print(String.valueOf(pixelPos.x));
        out.print("\t");
        out.print(String.valueOf(pixelPos.y));
        if (geoCoding != null) {
            final GeoPos geoPos = geoCoding.getGeoPos(pixelPos, null);
            out.print("\t");
            out.print(String.valueOf(geoPos.lon));
            out.print("\t");
            out.print(String.valueOf(geoPos.lat));
        }
        for (final Band band : bands) {
            out.print("\t");
            if (band.isFloatingPointType()) {
                band.readPixels(x, y, 1, 1, floatPixel, ProgressMonitor.NULL);
                out.print(floatPixel[0]);
            } else {
                band.readPixels(x, y, 1, 1, intPixel, ProgressMonitor.NULL);
                out.print(intPixel[0]);
            }
        }
        for (final TiePointGrid grid : tiePointGrids) {
            grid.readPixels(x, y, 1, 1, floatPixel, ProgressMonitor.NULL);
            out.print("\t");
            out.print(floatPixel[0]);
        }
        out.print("\n");
    }


    /**
     * Computes the total number of pixels within the specified ROI.
     *
     * @param raster   the raster data node
     * @param roiMaskImage the rendered image masking out the ROI
     * @return the total number of pixels in the ROI
     */
    private static long getNumROIPixels(final RasterDataNode raster, final RenderedImage roiMaskImage) {
        final int minTileX = roiMaskImage.getMinTileX();
        final int minTileY = roiMaskImage.getMinTileY();

        final int numXTiles = roiMaskImage.getNumXTiles();
        final int numYTiles = roiMaskImage.getNumYTiles();

        final int w = raster.getSceneRasterWidth();
        final int h = raster.getSceneRasterHeight();
        final Rectangle imageRect = new Rectangle(0, 0, w, h);

        long numROIPixels = 0;
        for (int tileX = minTileX; tileX < minTileX + numXTiles; ++tileX) {
            for (int tileY = minTileY; tileY < minTileY + numYTiles; ++tileY) {
                final Rectangle tileRectangle = new Rectangle(roiMaskImage.getTileGridXOffset() + tileX * roiMaskImage.getTileWidth(),
                                                              roiMaskImage.getTileGridYOffset() + tileY * roiMaskImage.getTileHeight(),
                                                              roiMaskImage.getTileWidth(), roiMaskImage.getTileHeight());

                final Rectangle r = imageRect.intersection(tileRectangle);
                if (!r.isEmpty()) {
                    Raster roiTile = roiMaskImage.getTile(tileX, tileY);
                    for (int y = r.y; y < r.y + r.height; y++) {
                        for (int x = r.x; x < r.x + r.width; x++) {
                            if (roiTile.getSample(x, y, 0) != 0) {
                                numROIPixels++;
                            }
                        }
                    }
                }
            }
        }
        return numROIPixels;
    }

    /**
     * Checks whether the command can be performed or not.
     *
     * @return <code>true</code> if so
     */
    private static boolean hasSelectedRasterMasks() {
        final ProductSceneView view = VisatApp.getApp().getSelectedProductSceneView();
        boolean enabled = false;

        if (view != null) {
            final RasterDataNode raster = view.getRaster();
            if (raster != null) {
                Product product = raster.getProduct();
                int numMasks = product.getMaskGroup().getNodeCount();
                enabled = numMasks > 0;
            }
        }
        return enabled;
    }

}
