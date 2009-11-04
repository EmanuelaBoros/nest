/*
 * $Id: SingleTargetProductDialog.java,v 1.2 2009-11-04 17:04:32 lveci Exp $
 *
 * Copyright (C) 2007 by Brockmann Consult (info@brockmann-consult.de)
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
package org.esa.beam.framework.gpf.ui;

import com.bc.ceres.core.ProgressMonitor;
import com.bc.ceres.core.SubProgressMonitor;
import com.bc.ceres.swing.progress.ProgressMonitorSwingWorker;
import org.esa.beam.framework.dataio.ProductIO;
import org.esa.beam.framework.datamodel.Product;
import org.esa.beam.framework.gpf.OperatorException;
import org.esa.beam.framework.gpf.operators.common.WriteOp;
import org.esa.beam.framework.ui.AppContext;
import org.esa.beam.framework.ui.BasicApp;
import org.esa.beam.framework.ui.ModelessDialog;
import org.esa.beam.framework.ui.SuppressibleOptionPane;
import org.esa.beam.util.SystemUtils;
import org.esa.beam.util.io.FileUtils;

import javax.swing.*;
import javax.media.jai.JAI;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.File;
import java.text.MessageFormat;
import java.util.concurrent.ExecutionException;
import java.util.Calendar;
import java.util.Date;
import java.awt.*;

/**
 * WARNING: This class belongs to a preliminary API and may change in future releases.
 *
 * @author Norman Fomferra
 * @author Marco Peters
 * @version $Revision: 1.2 $ $Date: 2009-11-04 17:04:32 $
 */
public abstract class SingleTargetProductDialog extends ModelessDialog {

    private TargetProductSelector targetProductSelector;
    private AppContext appContext;
    private JLabel statusLabel;

    public SingleTargetProductDialog(AppContext appContext, String title, String helpID) {
        this(appContext, title, ID_APPLY_CLOSE_HELP, helpID);
    }

    public SingleTargetProductDialog(AppContext appContext, String title, int buttonMask, String helpID) {
        super(appContext.getApplicationWindow(), title, buttonMask, helpID);
        this.appContext = appContext;
        targetProductSelector = new TargetProductSelector();
        String homeDirPath = SystemUtils.getUserHomeDir().getPath();
        String saveDir = appContext.getPreferences().getPropertyString(BasicApp.PROPERTY_KEY_APP_LAST_SAVE_DIR, homeDirPath);
        targetProductSelector.getModel().setProductDir(new File(saveDir));
        targetProductSelector.getOpenInAppCheckBox().setText("Open in " + appContext.getApplicationName());
        targetProductSelector.getModel().getValueContainer().addPropertyChangeListener(new PropertyChangeListener() {

            @Override
            public void propertyChange(PropertyChangeEvent evt) {
                if (evt.getPropertyName().equals("saveToFileSelected") ||
                        evt.getPropertyName().equals("openInAppSelected")) {
                    updateRunButton();
                }
            }
        });
        AbstractButton button = getButton(ID_APPLY);
        button.setText("Run");
        button.setMnemonic('R');
        updateRunButton();

        statusLabel = new JLabel("");
        statusLabel.setForeground(new Color(255,0,0));
        this.getJDialog().getContentPane().add(statusLabel, BorderLayout.NORTH);
    }

    private void updateRunButton() {
        AbstractButton button = getButton(ID_APPLY);
        boolean save = targetProductSelector.getModel().isSaveToFileSelected();
        boolean open = targetProductSelector.getModel().isOpenInAppSelected();
        String toolTipText = "";
        boolean enabled = true;
        if (save && open) {
            toolTipText = "Save target product and open it in " + getAppContext().getApplicationName();
        } else if (save) {
            toolTipText = "Save target product";
        } else if (open) {
            toolTipText = "Open target product in " + getAppContext().getApplicationName();
        } else {
            enabled = false;
        }
        button.setToolTipText(toolTipText);
        button.setEnabled(enabled);
    }

    public AppContext getAppContext() {
        return appContext;
    }

    public TargetProductSelector getTargetProductSelector() {
        return targetProductSelector;
    }

    @Override
    protected void onApply() {
        if (!canApply()) {
            return;
        }

        String productDir = targetProductSelector.getModel().getProductDir().getAbsolutePath();
        appContext.getPreferences().setPropertyString(BasicApp.PROPERTY_KEY_APP_LAST_SAVE_DIR, productDir);

        Product targetProduct = null;
        try {
            targetProduct = createTargetProduct();
            if (targetProduct == null) {
                throw new NullPointerException("Target product is null.");
            }
        } catch (Throwable t) {
            handleInitialisationError(t);
        }
        if (targetProduct == null) {
            return;
        }

        targetProduct.setName(targetProductSelector.getModel().getProductName());
        if (targetProductSelector.getModel().isSaveToFileSelected()) {
            targetProduct.setFileLocation(targetProductSelector.getModel().getProductFile());
            final ProgressMonitorSwingWorker worker = new ProductWriterSwingWorker(targetProduct);
            //worker.executeWithBlocking();
            worker.execute();
        } else if (targetProductSelector.getModel().isOpenInAppSelected()) {
            appContext.getProductManager().addProduct(targetProduct);
            showOpenInAppInfo();
        }
    }

    protected void handleInitialisationError(Throwable t) {
        String msg;
        if (t instanceof RuntimeException || t instanceof Error) {
            msg = MessageFormat.format("An internal error occurred during the target product initialisation.\n{0}",
                                       formatThrowable(t));
        } else {
            msg = MessageFormat.format("An error occurred during the target product initialisation.\n{0}",
                                       formatThrowable(t));
        }
        appContext.handleError(msg, t);
    }

    protected void handleProcessingError(Throwable t) {
        String msg;
        if (t instanceof RuntimeException || t instanceof Error) {
            msg = MessageFormat.format("An internal error occurred during the target product processing.\n{0}",
                                       formatThrowable(t));
        } else {
            msg = MessageFormat.format("An error occurred during processing the target product processing.\n{0}",
                                       formatThrowable(t));
        }
        appContext.handleError(msg, t);
    }

    private String formatThrowable(Throwable t) {
        return MessageFormat.format("Type: {0}\nMessage: {1}\n",
                                    t.getClass().getSimpleName(),
                                    t.getMessage());
    }

    private boolean canApply() {
        final String productName = targetProductSelector.getModel().getProductName();
        if (productName == null || productName.isEmpty()) {
            showErrorDialog("Please specify a target product name.");
            targetProductSelector.getProductNameTextField().requestFocus();
            return false;
        }

        if (targetProductSelector.getModel().isOpenInAppSelected()) {
            final Product existingProduct = appContext.getProductManager().getProduct(productName);
            if (existingProduct != null) {
                String message = MessageFormat.format(
                        "A product with the name ''{0}'' is already opened in {1}.\n\n" +
                                "Do you want to continue?",
                        productName, appContext.getApplicationName());
                final int answer = JOptionPane.showConfirmDialog(getJDialog(), message,
                                                                 getTitle(), JOptionPane.YES_NO_OPTION);
                if (answer != JOptionPane.YES_OPTION) {
                    return false;
                }
            }
        }
        if (targetProductSelector.getModel().isSaveToFileSelected()) {
            File productFile = targetProductSelector.getModel().getProductFile();
            if (productFile.exists()) {
                String message = MessageFormat.format(
                        "The specified output file\n\"{0}\"\n already exists.\n\n" +
                                "Do you want to overwrite the existing file?",
                        productFile.getPath());
                final int answer = JOptionPane.showConfirmDialog(getJDialog(), message,
                                                                 getTitle(), JOptionPane.YES_NO_OPTION);
                if (answer != JOptionPane.YES_OPTION) {
                    return false;
                }
            }
        }
        return true;
    }

    private void showSaveInfo(long saveTime) {
        File productFile = getTargetProductSelector().getModel().getProductFile();
        final String message = MessageFormat.format(
                "The target product has been successfully written to\n{0}\n" +
                        "Total time spend for processing: {2}",
                formatFile(productFile),
                formatDuration(saveTime));
        showSuppressibleInformationDialog(message, "saveInfo");
    }

    private void showOpenInAppInfo() {
        final String message = MessageFormat.format(
                "The target product has successfully been created and opened in {0}.\n\n" +
                        "Actual processing of source to target data will be performed only on demand,\n" +
                        "for example, if the target product is saved or an image view is opened.",
                appContext.getApplicationName());
        showSuppressibleInformationDialog(message, "openInAppInfo");
    }

    private void showSaveAndOpenInAppInfo(long saveTime) {
        File productFile = getTargetProductSelector().getModel().getProductFile();
        final String message = MessageFormat.format(
                "The target product has been successfully written to\n" +
                        "{0}\n" +
                        "and has been opened in {1}.\n" +
                        "Total time spend for processing: {2}\n",
                formatFile(productFile),
                appContext.getApplicationName(),
                formatDuration(saveTime));
        showSuppressibleInformationDialog(message, "saveAndOpenInAppInfo");
    }

    String formatFile(File file) {
        return FileUtils.getDisplayText(file, 54);
    }

    String formatDuration(long millis) {
        long seconds = millis / 1000;
        millis -= seconds * 1000;
        long minutes = seconds / 60;
        seconds -= minutes * 60;
        long hours = minutes / 60;
        minutes -= hours * 60;
        return String.format("%02d:%02d:%02d.%03d", hours, minutes, seconds, millis);
    }

    /**
     * Shows an information dialog on top of this dialog.
     * The shown dialog will contain a user option allowing to not show the message anymore.
     *
     * @param infoMessage  The message.
     * @param propertyName The (simple) property name used to store the user option in the application's preferences.
     */
    public void showSuppressibleInformationDialog(String infoMessage, String propertyName) {
        final SuppressibleOptionPane optionPane = new SuppressibleOptionPane(appContext.getPreferences());
        optionPane.showMessageDialog(getClass().getName() + "." + propertyName,
                                     getJDialog(),
                                     infoMessage,
                                     getTitle(),
                                     JOptionPane.INFORMATION_MESSAGE);
    }

    private class ProductWriterSwingWorker extends ProgressMonitorSwingWorker<Product, Object> {
        private final Product targetProduct;
        private long saveTime;
        private Date executeStartTime;

        private ProductWriterSwingWorker(Product targetProduct) {
            super(getJDialog(), "Writing Target Product");
            this.targetProduct = targetProduct;
        }

        @Override
        protected Product doInBackground(ProgressMonitor pm) throws Exception {
            final TargetProductSelectorModel model = getTargetProductSelector().getModel();
            pm.beginTask("Writing...", model.isOpenInAppSelected() ? 100 : 95);
            saveTime = 0L;
            Product product = null;
            try {
                JAI.getDefaultInstance().getTileCache().flush();
                System.gc();

                executeStartTime = Calendar.getInstance().getTime();
                long t0 = System.currentTimeMillis();
                WriteOp.writeProduct(targetProduct,
                                     model.getProductFile(),
                                     model.getFormatName(), SubProgressMonitor.create(pm, 95));
                saveTime = System.currentTimeMillis() - t0;
                if (model.isOpenInAppSelected()) {
                    product = ProductIO.readProduct(model.getProductFile(), null);
                    if (product == null) {
                        product = targetProduct; // todo - check - this cannot be ok!!! (nf)
                    }
                    pm.worked(5);
                }
            } finally {
                pm.done();
                if (product != targetProduct) {
                    targetProduct.dispose();
                }
            }
            return product;
        }

        @Override
        protected void done() {
            final TargetProductSelectorModel model = getTargetProductSelector().getModel();
            try {
                final Date now = Calendar.getInstance().getTime();
                final long diff = (now.getTime() - executeStartTime.getTime()) / 1000;
                if(diff > 120) {
                    final float minutes = diff / 60f;
                    statusLabel.setText("Processing completed in " + minutes + " minutes");
                } else {
                    statusLabel.setText("Processing completed in " + diff + " seconds");
                }

                final Product targetProduct = get();
                if (model.isOpenInAppSelected()) {
                    appContext.getProductManager().addProduct(targetProduct);
                    showSaveAndOpenInAppInfo(saveTime);
                } else {
                    showSaveInfo(saveTime);
                }
            } catch (InterruptedException e) {
                // ignore
            } catch (ExecutionException e) {
                handleProcessingError(e.getCause());
            } catch (Throwable t) {
                handleProcessingError(t);
            }
        }
    }

    /**
     * Creates the desired target product.
     * Usually, this method will be implemented by invoking one of the multiple {@link org.esa.beam.framework.gpf.GPF GPF}
     * {@code createProduct} methods.
     * <p/>
     * The method should throw a {@link OperatorException} in order to signal "nominal" processing errors,
     * other exeption types are treated as internal errors.
     *
     * @return The target product.
     * @throws Exception if an error occurs, an {@link OperatorException} is signaling "nominal" processing errors.
     */
    protected abstract Product createTargetProduct() throws Exception;
}