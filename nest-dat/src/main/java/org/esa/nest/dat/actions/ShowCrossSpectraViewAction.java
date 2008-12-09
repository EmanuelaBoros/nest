
package org.esa.nest.dat.actions;

import org.esa.beam.framework.datamodel.*;
import org.esa.beam.framework.ui.command.CommandEvent;
import org.esa.beam.framework.ui.command.ExecCommand;
import org.esa.beam.framework.ui.UIUtils;
import org.esa.beam.framework.ui.product.ProductSceneImage;
import org.esa.beam.framework.ui.product.ProductSceneView;
import org.esa.beam.visat.VisatApp;
import org.esa.beam.util.Debug;
import org.esa.nest.dat.views.polarview.PolarView;

import javax.swing.*;
import javax.swing.event.InternalFrameAdapter;
import javax.swing.event.InternalFrameEvent;

import com.bc.ceres.swing.progress.ProgressMonitorSwingWorker;
import com.bc.ceres.core.*;
import com.bc.ceres.core.ProgressMonitor;

/**
 * This action opens a polar wave view for the currently selected wave product.
 *
 */
public class ShowCrossSpectraViewAction extends ExecCommand {
    public static String ID = "showPolarWaveView";

    @Override
    public void actionPerformed(final CommandEvent event) {
        final VisatApp visatApp = VisatApp.getApp();
        openProductSceneView((RasterDataNode) visatApp.getSelectedProductNode());
    }

    public void openProductSceneView(final RasterDataNode selectedProductNode) {
        final VisatApp visatApp = VisatApp.getApp();
        visatApp.setStatusBarMessage("Creating polar view...");
        UIUtils.setRootFrameWaitCursor(visatApp.getMainFrame());
        final Product product = VisatApp.getApp().getSelectedProduct();

        final SwingWorker worker = new ProgressMonitorSwingWorker<ProductSceneImage, Object>(visatApp.getMainFrame(),
               visatApp.getAppName() + " - Creating image for '" + selectedProductNode.getName() + "'") {

                @Override
                protected ProductSceneImage doInBackground(ProgressMonitor pm) throws Exception {
                    try {
                        return createProductSceneImage(selectedProductNode, pm);
                    } finally {
                        if (pm.isCanceled()) {
                            selectedProductNode.unloadRasterData();
                        }
                    }
                }

                @Override
                public void done() {
                    UIUtils.setRootFrameDefaultCursor(visatApp.getMainFrame());
                    visatApp.clearStatusBarMessage();

                    final ProductSceneImage productSceneImage;
                    try {
                        productSceneImage = get();
                    } catch (OutOfMemoryError e) {
                        visatApp.showOutOfMemoryErrorDialog("The polar view could not be created.");
                        return;
                    } catch (Exception e) {
                        visatApp.handleUnknownException(e);
                        return;
                    }

                    PolarView view = new PolarView(product, productSceneImage);

                    final String title = createInternalFrameTitle(selectedProductNode);
                    final Icon icon = UIUtils.loadImageIcon("icons/RsBandAsSwath16.gif");
                    final JInternalFrame internalFrame = visatApp.createInternalFrame(title, icon, view, getHelpId());
                    final ProductNodeListenerAdapter pnl = new ProductNodeListenerAdapter() {
                        @Override
                        public void nodeChanged(final ProductNodeEvent event1) {
                            if (event1.getSourceNode() == selectedProductNode &&
                                    event1.getPropertyName().equalsIgnoreCase(ProductNode.PROPERTY_NAME_NAME)) {
                                internalFrame.setTitle(createInternalFrameTitle(selectedProductNode));
                            }
                        }
                    };
                    final Product product = selectedProductNode.getProduct();
                    internalFrame.addInternalFrameListener(new InternalFrameAdapter() {
                        @Override
                        public void internalFrameOpened(InternalFrameEvent event1) {
                            product.addProductNodeListener(pnl);
                        }

                        @Override
                        public void internalFrameClosed(InternalFrameEvent event11) {
                            product.removeProductNodeListener(pnl);
                        }
                    });

                    visatApp.updateState();
                }
            };
        visatApp.getExecutorService().submit(worker);
    }

    private String createInternalFrameTitle(final RasterDataNode raster) {
        return UIUtils.getUniqueFrameTitle(VisatApp.getApp().getAllInternalFrames(), raster.getDisplayName());
    }

    private ProductSceneImage createProductSceneImage(final RasterDataNode raster,
                                                     ProgressMonitor pm) {
        Debug.assertNotNull(raster);
        Debug.assertNotNull(pm);
        final VisatApp app = VisatApp.getApp();

        ProductSceneImage sceneImage = null;
        try {
            pm.beginTask("Creating polar view...", 1);
            final JInternalFrame[] frames = app.findInternalFrames(raster, 1);
            if (frames.length > 0) {
                final ProductSceneView view = (ProductSceneView) frames[0].getContentPane();
                sceneImage = new ProductSceneImage(raster, view);
            } else {
                sceneImage = new ProductSceneImage(raster, app.getPreferences(), SubProgressMonitor.create(pm, 1));
            }
        } finally {
            pm.done();
        }

        return sceneImage;
    }

    @Override
    public void updateState(final CommandEvent event) {
        Product product = VisatApp.getApp().getSelectedProduct();
        if(product != null) {
            final String productType = VisatApp.getApp().getSelectedProduct().getProductType();
            setEnabled(productType.startsWith("ASA_WV") &&
                VisatApp.getApp().getSelectedProductNode() instanceof RasterDataNode);
        } else
            setEnabled(false);
    }
}