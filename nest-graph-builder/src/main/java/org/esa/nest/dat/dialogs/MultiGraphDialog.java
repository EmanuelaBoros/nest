package org.esa.nest.dat.dialogs;

import com.bc.ceres.core.ProgressMonitor;
import com.bc.ceres.core.SubProgressMonitor;
import org.esa.beam.framework.dataio.ProductIO;
import org.esa.beam.framework.datamodel.Product;
import org.esa.beam.framework.gpf.graph.GraphException;
import org.esa.beam.framework.ui.AppContext;
import org.esa.beam.framework.ui.ModelessDialog;
import org.esa.nest.dat.plugins.graphbuilder.GraphExecuter;
import org.esa.nest.dat.plugins.graphbuilder.ProgressBarProgressMonitor;

import javax.media.jai.JAI;
import javax.swing.*;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Date;

/**
 *  Provides the dialog for excuting multiple graph from one user interface
 */
public abstract class MultiGraphDialog extends ModelessDialog {

    protected final AppContext appContext;
    protected final IOPanel ioPanel;
    protected final ArrayList<GraphExecuter> graphExecuterList = new ArrayList<GraphExecuter>(3);

    private final JPanel mainPanel;
    protected final JTabbedPane tabbedPane;
    private final JLabel statusLabel;
    private final JPanel progressPanel;
    private final JProgressBar progressBar;
    private ProgressBarProgressMonitor progBarMonitor = null;

    private boolean isProcessing = false;

    protected static final String TMP_FILENAME = "tmp_intermediate";

    public MultiGraphDialog(final AppContext theAppContext, final String title, final String helpID,
                            final boolean useSourceSelector) {
        super(theAppContext.getApplicationWindow(), title, ID_APPLY_CLOSE_HELP, helpID);
        appContext = theAppContext;

        mainPanel = new JPanel(new BorderLayout(4, 4));

        tabbedPane = new JTabbedPane();
        tabbedPane.addChangeListener(new ChangeListener() {

            public void stateChanged(final ChangeEvent e) {
                ValidateAllNodes();
            }
        });
        mainPanel.add(tabbedPane, BorderLayout.CENTER);

        ioPanel = new IOPanel(appContext, tabbedPane, useSourceSelector);

        // status
        statusLabel = new JLabel("");
        statusLabel.setForeground(new Color(255,0,0));
        mainPanel.add(statusLabel, BorderLayout.NORTH);

        // progress Bar
        progressBar = new JProgressBar();
        progressBar.setName(getClass().getName() + "progressBar");
        progressBar.setStringPainted(true);
        progressPanel = new JPanel();
        progressPanel.setLayout(new BorderLayout(2,2));
        progressPanel.add(progressBar, BorderLayout.CENTER);
        final JButton progressCancelBtn = new JButton("Cancel");
        progressCancelBtn.addActionListener(new ActionListener() {

            public void actionPerformed(final ActionEvent e) {
                CancelProcessing();
            }
        });
        progressPanel.add(progressCancelBtn, BorderLayout.EAST);
        progressPanel.setVisible(false);
        mainPanel.add(progressPanel, BorderLayout.SOUTH);

        getButton(ID_APPLY).setText("Run");

        super.getJDialog().setMinimumSize(new Dimension(400, 300));
    }

    @Override
    public int show() {
        ioPanel.initProducts();
        setContent(mainPanel);
        initGraphs();
        return super.show();
    }

    @Override
    public void hide() {
        ioPanel.releaseProducts();
        super.hide();
    }

    @Override
    protected void onApply() {

        if(isProcessing) return;

        ioPanel.onApply();

        try {
            DoProcessing();
        } catch(Exception e) {
            statusLabel.setText(e.getMessage());
        }
    }

    @Override
    protected void onClose() {
        CancelProcessing();
        
        super.onClose();
    }

    void initGraphs() {
        try {
            deleteGraphs();
            createGraphs();
        } catch(Exception e) {
            statusLabel.setText(e.getMessage());
        }
    }

    /**
     * Validates the input and then call the GPF to execute the graph
     * @throws GraphException on assignParameters
     */
    private void DoProcessing() throws GraphException {

        if(ValidateAllNodes()) {

            JAI.getDefaultInstance().getTileCache().flush();
            System.gc();

            progressBar.setValue(0);
            progBarMonitor = new ProgressBarProgressMonitor(progressBar, null, progressPanel);

            final SwingWorker processThread = new ProcessThread(progBarMonitor);
            processThread.execute();

        } else {
            showErrorDialog(statusLabel.getText());
        }
    }

    private void CancelProcessing() {
        if(progBarMonitor != null)
            progBarMonitor.setCanceled(true);
    }

    private void deleteGraphs() {
        for(GraphExecuter gex : graphExecuterList) {
            gex.ClearGraph();
        }
        graphExecuterList.clear();
    }

    /**
     * Loads a new graph from a file
     * @param executer the GraphExcecuter
     * @param file the graph file to load
     */
    public void LoadGraph(final GraphExecuter executer, final File file) {
        try {
            executer.loadGraph(file, true);

        } catch(GraphException e) {
            showErrorDialog(e.getMessage());
        }
    }

    protected abstract void createGraphs() throws GraphException;

    protected abstract void assignParameters() throws GraphException;

    protected abstract void cleanUpTempFiles();

    private boolean ValidateAllNodes() {
        if(isProcessing) return false;
        if(ioPanel == null || graphExecuterList.isEmpty())
            return false;

        boolean result;
        statusLabel.setText("");
        try { 
            assignParameters();
            // first graph must pass
            result = graphExecuterList.get(0).InitGraph();

        } catch(GraphException e) {
            statusLabel.setText(e.getMessage());
            result = false;
        }
        return result;
    }

    private void openTargetProducts(final ArrayList<File> fileList) {
        if(!fileList.isEmpty()) {
            for(File file : fileList) {
                try {

                    final Product product = ProductIO.readProduct(file, null);
                    if (product != null) {
                        appContext.getProductManager().addProduct(product);
                    }
                } catch(IOException e) {
                    showErrorDialog(e.getMessage());
                }
            }
        }
    }

    protected IOPanel getIOPanel() {
        return ioPanel;
    }

    public void setTargetProductNameSuffix(final String suffix) {
        ioPanel.setTargetProductNameSuffix(suffix);
    }

    /////

    private class ProcessThread extends SwingWorker<Boolean, Object> {

        private final ProgressMonitor pm;
        private Date executeStartTime = null;

        public ProcessThread(final ProgressMonitor pm) {
            this.pm = pm;
        }

        @Override
        protected Boolean doInBackground() throws Exception {

            pm.beginTask("Processing Graph...", 100*graphExecuterList.size());
            try {
                executeStartTime = Calendar.getInstance().getTime();
                isProcessing = true;

                for(GraphExecuter graphEx : graphExecuterList) {
                    graphEx.InitGraph();

                    graphEx.executeGraph(new SubProgressMonitor(pm, 100));
                    graphEx.disposeGraphContext();
                }

            } catch(Exception e) {
                System.out.print(e.getMessage());
            } finally {
                isProcessing = false;
                pm.done();
            }
            return true;
        }

        @Override
        public void done() {
            final Date now = Calendar.getInstance().getTime();
            final long diff = (now.getTime() - executeStartTime.getTime()) / 1000;
            if(diff > 120) {
                final float minutes = diff / 60f;
                statusLabel.setText("Processing completed in " + minutes + " minutes");
            } else {
                statusLabel.setText("Processing completed in " + diff + " seconds");
            }

            if(ioPanel.isOpenInAppSelected()) {
                final GraphExecuter graphEx = graphExecuterList.get(graphExecuterList.size()-1);
                openTargetProducts(graphEx.getProductsToOpenInDAT());
            }

            cleanUpTempFiles();
        }

    }

}