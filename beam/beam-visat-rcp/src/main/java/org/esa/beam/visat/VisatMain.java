/*
 * $Id: VisatMain.java,v 1.3 2009-11-04 17:04:32 lveci Exp $
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
package org.esa.beam.visat;

import com.bc.ceres.core.ProgressMonitor;
import com.bc.ceres.core.runtime.RuntimeRunnable;
import com.jidesoft.utils.Lm;
import com.jidesoft.utils.SystemInfo;
import org.esa.beam.framework.dataio.ProductIO;
import org.esa.beam.framework.datamodel.Product;
import org.esa.beam.framework.ui.BasicApp;
import org.esa.beam.framework.ui.UIUtils;
import org.esa.beam.framework.ui.application.ApplicationDescriptor;
import org.esa.beam.util.Debug;
import org.esa.beam.visat.actions.session.OpenSessionAction;
import org.esa.beam.BeamUiActivator;

import javax.swing.JOptionPane;
import javax.swing.SwingUtilities;
import javax.media.jai.JAI;
import javax.media.jai.util.ImagingListener;
import java.io.File;
import java.io.IOException;
import java.text.MessageFormat;
import java.util.ArrayList;
import java.util.Locale;

/**
 * The startup class for VISAT. It provides the <code>main</code> method for the application.
 * <p/>
 * <p>The VISAT application accepts the following command line options: <ld> <li> <code>-d</code> or
 * <code>--debug</code> sets VISAT into debugging mode <li> <code>-l <i>file</i></code> or <code>--logfile
 * <i>file</i></code> sets the logfile for VISAT to <i>file</i> </ld>
 *
 * @author Norman Fomferra
 * @version $Revision: 1.3 $ $Date: 2009-11-04 17:04:32 $
 */
public class VisatMain implements RuntimeRunnable {
    /**
     * Entry point for the VISAT application called by the Ceres runtime.
     *
     * @param argument        a {@code String[]} containing the command line arguments
     * @param progressMonitor a progress monitor
     * @throws Exception if an error occurs
     */
    @Override
    public void run(Object argument, ProgressMonitor progressMonitor) throws Exception {

        String[] args = new String[0];
        if (argument instanceof String[]) {
            args = (String[]) argument;
        }

        Locale.setDefault(Locale.UK); // Force usage of British English locale

        Lm.verifyLicense("Brockmann Consult", "BEAM", "lCzfhklpZ9ryjomwWxfdupxIcuIoCxg2");
        if (SystemInfo.isMacOSX()) {
            if (System.getProperty("com.apple.macos.useScreenMenuBar") == null) {
                System.setProperty("com.apple.macos.useScreenMenuBar", "true");
            }
            if (System.getProperty("apple.laf.useScreenMenuBar") == null) {
                System.setProperty("apple.laf.useScreenMenuBar", "true");
            }
            if (System.getProperty("apple.awt.brushMetalLook") == null) {
                System.setProperty("apple.awt.brushMetalLook", "true");
            }
        }

        final ApplicationDescriptor applicationDescriptor = BeamUiActivator.getInstance().getApplicationDescriptor();
        if (applicationDescriptor == null) {
             throw new IllegalStateException(String.format("Application descriptor not found for applicationId='%s'.",
                                                           BeamUiActivator.getInstance().getApplicationId()));
        }

        boolean debugEnabled = false;
        ArrayList<String> productFilepathList = new ArrayList<String>();
        String sessionFile = null;
        for (String arg : args) {
            if (arg.startsWith("-")) {
                if (arg.equals("-d") || arg.equals("--debug")) {
                    debugEnabled = true;
                } else {
                    System.err.printf("%s error: illegal option '" + arg + "'", applicationDescriptor.getDisplayName());
                    return;
                }
            } else if (arg.endsWith(OpenSessionAction.SESSION_FILE_FILTER.getDefaultExtension())) {
                sessionFile = arg;
            } else {
                productFilepathList.add(arg);
            }
        }

        Debug.setEnabled(debugEnabled);
        if (debugEnabled) {
            JAI.getDefaultInstance().setImagingListener(new ImagingListener() {
                @Override
                public boolean errorOccurred(String message, Throwable thrown, Object where, boolean isRetryable) throws RuntimeException {
                    Debug.trace("JAI Error: " + message);
                    Debug.trace(thrown);
                    return false;
                }
            });
        }
        JAI.getDefaultInstance().getTileScheduler().setParallelism(Runtime.getRuntime().availableProcessors());

        final VisatApp app = createApplication(applicationDescriptor);
        app.startUp(progressMonitor);
        openSession(app, sessionFile);
        openProducts(app, productFilepathList);
    }

    protected VisatApp createApplication(ApplicationDescriptor applicationDescriptor) {
        return new VisatApp(applicationDescriptor);
    }
    
    private void openSession(VisatApp app, String sessionFile) {
        if (sessionFile != null && !(sessionFile.trim().isEmpty())) {
            final OpenSessionAction action = (OpenSessionAction) app.getCommandManager().getCommand(OpenSessionAction.ID);
            action.openSession(app, new File(sessionFile));
        }
    }

    private static void openProducts(VisatApp app, ArrayList<String> productFilepathList) {
        for (String productFilepath : productFilepathList) {
            openProduct(app, productFilepath);
        }
    }

    private static void openProduct(final VisatApp app, final String productFilepath) {
        SwingUtilities.invokeLater(new Runnable() {
            @Override
            public void run() {
                UIUtils.setRootFrameWaitCursor(app.getMainFrame());
                try {
                    openProductImpl(app, productFilepath);
                } finally {
                    UIUtils.setRootFrameDefaultCursor(app.getMainFrame());
                }
            }
        });
    }

    private static void openProductImpl(VisatApp app, final String productFilepath) {
        final File productFile = new File(productFilepath);
        final Product product;
        try {
            product = ProductIO.readProduct(productFile, null);
            if (product == null) {
                final MessageFormat mf = new MessageFormat("No reader found for data product\n''{0}''."); /*I18N*/
                final Object[] args = new Object[]{productFile.getPath()};
                showError(app, mf.format(args));
                return;
            }
        } catch (IOException e) {
            final MessageFormat mf = new MessageFormat("I/O error while opening file\n{0}:\n{1}"); /*I18N*/
            final Object[] args = new Object[]{productFile.getPath(), e.getMessage()};
            showError(app, mf.format(args));
            return;
        }
        app.addProduct(product);
    }

    private static void showError(BasicApp app, final String message) {
        JOptionPane.showMessageDialog(null,
                                      message,
                                      app.getAppName(),
                                      JOptionPane.ERROR_MESSAGE);
    }
}
