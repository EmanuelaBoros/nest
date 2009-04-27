package org.esa.beam.visat.toolviews.stat;

import com.bc.ceres.core.ProgressMonitor;
import com.bc.ceres.core.SubProgressMonitor;
import com.bc.ceres.swing.TableLayout;
import com.bc.ceres.swing.progress.ProgressMonitorSwingWorker;
import org.esa.beam.framework.datamodel.Band;
import org.esa.beam.framework.datamodel.Product;
import org.esa.beam.framework.datamodel.RasterDataNode;
import org.esa.beam.framework.datamodel.Stx;
import org.esa.beam.framework.datamodel.TiePointGrid;
import org.esa.beam.framework.param.ParamChangeEvent;
import org.esa.beam.framework.param.ParamChangeListener;
import org.esa.beam.framework.param.ParamGroup;
import org.esa.beam.framework.param.Parameter;
import org.esa.beam.framework.param.editors.ComboBoxEditor;
import org.esa.beam.framework.ui.GridBagUtils;
import org.esa.beam.framework.ui.application.ToolView;
import org.esa.beam.util.Debug;
import org.esa.beam.util.ProductUtils;
import org.esa.beam.util.math.MathUtils;
import org.esa.beam.util.math.Range;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.title.TextTitle;
import org.jfree.chart.title.Title;
import org.jfree.ui.RectangleInsets;

import javax.swing.BorderFactory;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.GridBagConstraints;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.awt.image.DataBufferByte;
import java.awt.image.IndexColorModel;
import java.awt.image.RenderedImage;
import java.io.IOException;
import java.text.MessageFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.CancellationException;
import java.util.concurrent.ExecutionException;


/**
 * The scatter plot pane within the statistcs window.
 *
 * @author Marco Peters
 */
class ScatterPlotPanel extends PagePanel {

    private static final String NO_DATA_MESSAGE = "No scatter plot computed yet.\n" +
                                                  "TIP: To zoom within the chart draw a rectangle\n" +
                                                  "with the mouse or use the context menu.";  /*I18N*/
    private static final String CHART_TITLE = "Scatter Plot";
    private static final String TITLE_PREFIX = CHART_TITLE;

    private static final int X_VAR = 0;
    private static final int Y_VAR = 1;

    private ParamGroup paramGroup;

    private static Parameter[] rasterNameParams = new Parameter[2];
    private static Parameter[] autoMinMaxParams = new Parameter[2];
    private static Parameter[] minParams = new Parameter[2];
    private static Parameter[] maxParams = new Parameter[2];

    private ComputePanel computePanel;
    private ChartPanel scatterPlotDisplay;
    private boolean adjustingAutoMinMax;
    private XYImagePlot plot;

    ScatterPlotPanel(ToolView parentDialog, String helpId) {
        super(parentDialog, helpId);
    }

    @Override
    public String getTitle() {
        return getTitlePrefix();
    }

    @Override
    protected String getTitlePrefix() {
        return TITLE_PREFIX;
    }

    @Override
    protected void initContent() {
        initParameters();
        createUI();
        updateContent();
    }

    @Override
    protected void updateContent() {
        if (scatterPlotDisplay != null) {
            final String[] availableBands = createAvailableBandList();
            updateParameters(X_VAR, availableBands);
            updateParameters(Y_VAR, availableBands);
            setChartTitle();
        }
    }

    private void setChartTitle() {
        final JFreeChart chart = scatterPlotDisplay.getChart();
        final List<Title> subtitles = new ArrayList<Title>(7);
        subtitles.add(new TextTitle(MessageFormat.format("{0}, {1}",
                                                         rasterNameParams[X_VAR].getValueAsText(),
                                                         rasterNameParams[Y_VAR].getValueAsText())));
        chart.setSubtitles(subtitles);
    }

    private void updateParameters(int varIndex, String[] availableBands) {

        final RasterDataNode raster = getRaster();
        String rasterName = null;
        if (raster != null) {
            rasterName = raster.getName();
        } else if (availableBands.length > 0) {
            rasterName = availableBands[0];
        }
        if (rasterName != null) {
            rasterNameParams[varIndex].getProperties().setValueSet(availableBands);
            rasterNameParams[varIndex].setValue(rasterName, null);
        } else {
            rasterNameParams[varIndex].getProperties().setValueSet(new String[0]);
        }

        if ((Boolean) autoMinMaxParams[varIndex].getValue()) {
            minParams[varIndex].setDefaultValue();
            maxParams[varIndex].setDefaultValue();
        }
    }

    private void initParameters() {
        paramGroup = new ParamGroup();

        final String[] availableBands = createAvailableBandList();
        initParameters(X_VAR, availableBands);
        initParameters(Y_VAR, availableBands);

        paramGroup.addParamChangeListener(new ParamChangeListener() {

            public void parameterValueChanged(ParamChangeEvent event) {
                updateUIState();
            }
        });
    }

    private void initParameters(int varIndex, String[] availableBands) {

        final String paramPrefix = "var" + varIndex + ".";

        final RasterDataNode raster = getRaster();
        final String rasterName;
        if (raster != null) {
            rasterName = raster.getName();
        } else if (availableBands.length > 0) {
            rasterName = availableBands[0];
        } else {
            rasterName = "";
        }
        rasterNameParams[varIndex] = new Parameter(paramPrefix + "rasterName", rasterName);
        rasterNameParams[varIndex].getProperties().setValueSet(availableBands);
        rasterNameParams[varIndex].getProperties().setValueSetBound(true);
        rasterNameParams[varIndex].getProperties().setDescription("Band name"); /*I18N*/
        rasterNameParams[varIndex].getProperties().setEditorClass(ComboBoxEditor.class);
        paramGroup.addParameter(rasterNameParams[varIndex]);

        autoMinMaxParams[varIndex] = new Parameter(paramPrefix + "autoMinMax", Boolean.TRUE);
        autoMinMaxParams[varIndex].getProperties().setLabel("Auto min/max");
        autoMinMaxParams[varIndex].getProperties().setDescription("Automatically detect min/max");  /*I18N*/
        paramGroup.addParameter(autoMinMaxParams[varIndex]);

        minParams[varIndex] = new Parameter(paramPrefix + "min", 0.0);
        minParams[varIndex].getProperties().setLabel("Min:");
        minParams[varIndex].getProperties().setDescription("Minimum display value");    /*I18N*/
        minParams[varIndex].getProperties().setNumCols(7);
        paramGroup.addParameter(minParams[varIndex]);

        maxParams[varIndex] = new Parameter(paramPrefix + "max", 100.0);
        maxParams[varIndex].getProperties().setLabel("Max:");
        maxParams[varIndex].getProperties().setDescription("Maximum display value");    /*I18N*/
        maxParams[varIndex].getProperties().setNumCols(7);
        paramGroup.addParameter(maxParams[varIndex]);
    }

    private void createUI() {
        final ActionListener actionAll = new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                computeScatterPlot(false);
            }
        };
        final ActionListener actionROI = new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                computeScatterPlot(true);
            }
        };
        computePanel = ComputePanel.createComputePane(actionAll, actionROI, getRaster());

        plot = new XYImagePlot();

        plot.setAxisOffset(new RectangleInsets(5, 5, 5, 5));
        plot.setNoDataMessage(NO_DATA_MESSAGE);
        plot.getRenderer().setBaseToolTipGenerator(new XYPlotToolTipGenerator());

//        NumberAxis scaleAxis = new NumberAxis("Scale");
//        scaleAxis.setTickLabelFont(new Font("Dialog", Font.PLAIN, 7));
//        PaintScaleLegend legend = new PaintScaleLegend(new GrayPaintScale(), scaleAxis);
//        legend.setAxisLocation(AxisLocation.BOTTOM_OR_LEFT);
//        legend.setAxisOffset(5.0);
//        legend.setMargin(new RectangleInsets(5, 5, 5, 5));
//        legend.setPadding(new RectangleInsets(10, 10, 10, 10));
//        legend.setStripWidth(10);
//        legend.setPosition(RectangleEdge.RIGHT);

        JFreeChart chart = new JFreeChart(CHART_TITLE, plot);
        chart.removeLegend();
//         chart.addSubtitle(legend);

        scatterPlotDisplay = new ChartPanel(chart);
        scatterPlotDisplay.getPopupMenu().add(createCopyDataToClipboardMenuItem());

        final TableLayout rightPanelLayout = new TableLayout(1);
        final JPanel rightPanel = new JPanel(rightPanelLayout);
        rightPanelLayout.setTableFill(TableLayout.Fill.HORIZONTAL);
        rightPanelLayout.setRowWeightY(3, 1.0);
        rightPanelLayout.setCellFill(5, 1, TableLayout.Fill.NONE);
        rightPanelLayout.setCellAnchor(5, 1, TableLayout.Anchor.EAST);
        rightPanel.add(computePanel);
        rightPanel.add(createOptionsPane());
        rightPanel.add(createChartButtonPanel(scatterPlotDisplay));
        rightPanel.add(new JPanel());   // filler
        final JPanel helpPanel = new JPanel(new BorderLayout());
        helpPanel.add(getHelpButton(), BorderLayout.EAST);
        rightPanel.add(helpPanel);

        add(scatterPlotDisplay, BorderLayout.CENTER);
        add(rightPanel, BorderLayout.EAST);

        updateUIState();
    }

    private RasterDataNode getRaster(int varIndex) {
        final Product product = getProduct();
        if (product == null) {
            return null;
        }
        final String rasterName = rasterNameParams[varIndex].getValue().toString();
        RasterDataNode raster = product.getRasterDataNode(rasterName);
        if (raster == null) {
            // todo - "if the product doesn't have the raster, take some other?" - this is stupid (nf - 18.12.2007)
            if (getRaster() != null && rasterName.equalsIgnoreCase(getRaster().getName())) {
                raster = getRaster();
            }
        }
        Debug.assertTrue(raster != null);
        return raster;
    }

    private void updateUIState() {
        updateComputePane();
        updateUIState(X_VAR);
        updateUIState(Y_VAR);
        setChartTitle();
    }

    private void updateComputePane() {
        computePanel.setRaster(getRaster(X_VAR));
    }

    private void updateUIState(int varIndex) {
        final double min = ((Number) minParams[varIndex].getValue()).doubleValue();
        final double max = ((Number) maxParams[varIndex].getValue()).doubleValue();
        if (!adjustingAutoMinMax && min > max) {
            minParams[varIndex].setValue(max, null);
            maxParams[varIndex].setValue(min, null);
        }
        final boolean autoMinMaxEnabled = (Boolean) autoMinMaxParams[varIndex].getValue();
        minParams[varIndex].setUIEnabled(!autoMinMaxEnabled);
        maxParams[varIndex].setUIEnabled(!autoMinMaxEnabled);
    }


    private JPanel createOptionsPane() {
        final JPanel optionsPane = GridBagUtils.createPanel();
        final GridBagConstraints gbc = GridBagUtils.createConstraints("anchor=NORTHWEST,fill=BOTH");

        GridBagUtils.setAttributes(gbc, "gridy=1,weightx=1");
        GridBagUtils.addToPanel(optionsPane, createOptionsPane(X_VAR), gbc, "gridy=0,insets.top=0");
        GridBagUtils.addToPanel(optionsPane, createOptionsPane(Y_VAR), gbc, "gridy=1,insets.top=7");

        return optionsPane;
    }

    private JPanel createOptionsPane(int varIndex) {

        final JPanel optionsPane = GridBagUtils.createPanel();
        final GridBagConstraints gbc = GridBagUtils.createConstraints("anchor=WEST,fill=HORIZONTAL");

        GridBagUtils.setAttributes(gbc, "gridwidth=2,gridy=0,insets.top=0");
        GridBagUtils.addToPanel(optionsPane, rasterNameParams[varIndex].getEditor().getComponent(), gbc,
                                "gridx=0,weightx=1");

        GridBagUtils.setAttributes(gbc, "gridwidth=2,gridy=1,insets.top=4");
        GridBagUtils.addToPanel(optionsPane, autoMinMaxParams[varIndex].getEditor().getComponent(), gbc,
                                "gridx=0,weightx=1");

        GridBagUtils.setAttributes(gbc, "gridwidth=1,gridy=2,insets.top=2");
        GridBagUtils.addToPanel(optionsPane, minParams[varIndex].getEditor().getLabelComponent(), gbc,
                                "gridx=0,weightx=0");
        GridBagUtils.addToPanel(optionsPane, minParams[varIndex].getEditor().getComponent(), gbc, "gridx=1,weightx=1");

        GridBagUtils.setAttributes(gbc, "gridwidth=1,gridy=3,insets.top=2");
        GridBagUtils.addToPanel(optionsPane, maxParams[varIndex].getEditor().getLabelComponent(), gbc,
                                "gridx=0,weightx=0");
        GridBagUtils.addToPanel(optionsPane, maxParams[varIndex].getEditor().getComponent(), gbc, "gridx=1,weightx=1");

        optionsPane.setBorder(BorderFactory.createTitledBorder((varIndex == 0 ? "X" : "Y") + "-Band"));

        return optionsPane;
    }

    private void computeScatterPlot(boolean useROI) {
        try {
            computeScatterPlotImpl(useROI);
        } catch (IOException e) {
            JOptionPane.showMessageDialog(this,
                                          "An I/O error occured:\n" + e.getMessage(), /*I18N*/
                                          CHART_TITLE, /*I18N*/
                                          JOptionPane.ERROR_MESSAGE);
        }
    }

    private static class ScatterPlot {

        private final BufferedImage image;
        private final Range rangeX;
        private final Range rangeY;

        private ScatterPlot(BufferedImage image, Range rangeX, Range rangeY) {
            this.image = image;
            this.rangeX = rangeX;
            this.rangeY = rangeY;
        }

        public BufferedImage getImage() {
            return image;
        }

        public Range getRangeX() {
            return rangeX;
        }

        public Range getRangeY() {
            return rangeY;
        }
    }

    private void computeScatterPlotImpl(boolean useROI) throws IOException {

        final RasterDataNode rasterX = getRaster(X_VAR);
        final RasterDataNode rasterY = getRaster(Y_VAR);

        if (rasterX == null || rasterY == null) {
            return;
        }

        final RenderedImage roiImage;
        if (useROI) {
            roiImage = getRoiImage(rasterX);
        } else {
            roiImage = null;
        }
        ProgressMonitorSwingWorker<ScatterPlot, Object> swingWorker = new ProgressMonitorSwingWorker<ScatterPlot, Object>(
                this, "Computing scatter plot") {

            @Override
            protected ScatterPlot doInBackground(ProgressMonitor pm) throws Exception {
                pm.beginTask("Computing scatter plot...", 100);
                try {
                    final Range rangeX = getRange(X_VAR, rasterX, roiImage, SubProgressMonitor.create(pm, 15));
                    final Range rangeY = getRange(Y_VAR, rasterY, roiImage, SubProgressMonitor.create(pm, 15));
                    final BufferedImage image = ProductUtils.createScatterPlotImage(rasterX,
                                                                                    (float) rangeX.getMin(),
                                                                                    (float) rangeX.getMax(),
                                                                                    rasterY,
                                                                                    (float) rangeY.getMin(),
                                                                                    (float) rangeY.getMax(),
                                                                                    roiImage,
                                                                                    512,
                                                                                    512,
                                                                                    new Color(255, 255, 255, 0),
                                                                                    null,
                                                                                    SubProgressMonitor.create(pm, 70));
                    return new ScatterPlot(image, rangeX, rangeY);
                } finally {
                    pm.done();
                }
            }

            @Override
            public void done() {
                try {
                    final ScatterPlot scatterPlot = get();
                    double minX = scatterPlot.getRangeX().getMin();
                    double maxX = scatterPlot.getRangeX().getMax();
                    double minY = scatterPlot.getRangeY().getMin();
                    double maxY = scatterPlot.getRangeY().getMax();
                    if (minX > maxX || minY > maxY) {
                        JOptionPane.showMessageDialog(getParentDialogContentPane(),
                                                      "Failed to compute scatter plot.\n" +
                                                      "No Pixels considered..",
                                                      /*I18N*/
                                                      CHART_TITLE, /*I18N*/
                                                      JOptionPane.ERROR_MESSAGE);
                        plot.setDataset(null);
                        return;

                    }

                    if (MathUtils.equalValues(minX, maxX, 1.0e-4)) {
                        minX = Math.floor(minX);
                        maxX = Math.ceil(maxX);
                    }
                    if (MathUtils.equalValues(maxY, maxY, 1.0e-4)) {
                        minY = Math.floor(minY);
                        maxY = Math.ceil(maxY);
                    }
                    plot.setImage(scatterPlot.getImage());
                    plot.setImageDataBounds(new Rectangle2D.Double(minX, minY, maxX - minX, maxY - minY));
                    setAutoRange(X_VAR, new Range(minX, maxX));
                    setAutoRange(Y_VAR, new Range(minY, maxY));
                    plot.getDomainAxis().setLabel(getAxisLabel(getRaster(X_VAR)));
                    plot.getRangeAxis().setLabel(getAxisLabel(getRaster(Y_VAR)));
                } catch (InterruptedException e) {
                    e.printStackTrace();
                    JOptionPane.showMessageDialog(getParentDialogContentPane(),
                                                  "Failed to compute scatter plot.\n" +
                                                  "Calculation canceled.",
                                                  /*I18N*/
                                                  CHART_TITLE, /*I18N*/
                                                  JOptionPane.ERROR_MESSAGE);
                } catch (CancellationException e) {
                    e.printStackTrace();
                    JOptionPane.showMessageDialog(getParentDialogContentPane(),
                                                  "Failed to compute scatter plot.\n" +
                                                  "Calculation canceled.",
                                                  /*I18N*/
                                                  CHART_TITLE, /*I18N*/
                                                  JOptionPane.ERROR_MESSAGE);
                } catch (ExecutionException e) {
                    e.printStackTrace();
                    JOptionPane.showMessageDialog(getParentDialogContentPane(),
                                                  "Failed to compute scatter plot.\n" +
                                                  "An error occured:\n" +
                                                  e.getCause().getMessage(),
                                                  CHART_TITLE, /*I18N*/
                                                  JOptionPane.ERROR_MESSAGE);
                }
            }
        };
        swingWorker.execute();
    }

    private String getAxisLabel(RasterDataNode raster) {
        final String unit = raster.getUnit();
        if (unit != null) {
            return raster.getName() + " (" + unit + ")";
        }
        return raster.getName();
    }

    private Range getRange(int varIndex, RasterDataNode raster, RenderedImage roiImage, ProgressMonitor pm) throws
                                                                                                            IOException {
        final boolean autoMinMax = (Boolean) autoMinMaxParams[varIndex].getValue();
        if (autoMinMax) {
            Stx stx;
            if (roiImage == null) {
                stx = raster.getStx(false, pm);
            } else {
                stx = Stx.create(raster, roiImage, pm);
            }
            return new Range(raster.scale(stx.getMin()), raster.scale(stx.getMax()));
        } else {
            return new Range((Double) minParams[varIndex].getValue(), (Double) maxParams[varIndex].getValue());
        }
    }

    private void setAutoRange(int varIndex, Range range) {
        final boolean autoMinMax = (Boolean) autoMinMaxParams[varIndex].getValue();
        if (autoMinMax) {
            adjustingAutoMinMax = true;
            minParams[varIndex].setValueAsText(String.format("%7.2f", range.getMin()), null);
            maxParams[varIndex].setValueAsText(String.format("%7.2f", range.getMax()), null);
            adjustingAutoMinMax = false;
        }
    }

    private String[] createAvailableBandList() {
        final Product product = getProduct();
        final List<String> availableBandList = new ArrayList<String>(17);
        if (product != null) {
            for (int i = 0; i < product.getNumBands(); i++) {
                final Band band = product.getBandAt(i);
                availableBandList.add(band.getName());
            }
            for (int i = 0; i < product.getNumTiePointGrids(); i++) {
                final TiePointGrid grid = product.getTiePointGridAt(i);
                availableBandList.add(grid.getName());
            }
        }
        // if raster is only binded to the product and does not belong to it
        final RasterDataNode raster = getRaster();
        if (raster != null && raster.getProduct() == product) {
            final String rasterName = raster.getName();
            if (!availableBandList.contains(rasterName)) {
                availableBandList.add(rasterName);
            }
        }

        return availableBandList.toArray(new String[availableBandList.size()]);
    }

    @Override
    protected boolean checkDataToClipboardCopy() {
        final int warnLimit = 2000;
        final int excelLimit = 65536;
        final int numNonEmptyBins = getNumNonEmptyBins();
        if (numNonEmptyBins > warnLimit) {
            String excelNote = "";
            if (numNonEmptyBins > excelLimit - 100) {
                excelNote = "Note that e.g., Microsoft® Excel 2002 only supports a total of "
                            + excelLimit + " rows in a sheet.\n";   /*I18N*/
            }
            final String message = MessageFormat.format(
                    "This scatter plot diagram contains {0} non-empty bins.\n" +
                    "For each bin, a text data row containing an x, y and z value will be created.\n" +
                    "{1}\nPress ''Yes'' if you really want to copy this amount of data to the system clipboard.\n",
                    numNonEmptyBins, excelNote);
            final int status = JOptionPane.showConfirmDialog(this,
                                                             message, /*I18N*/
                                                             "Copy Data to Clipboard", /*I18N*/
                                                             JOptionPane.YES_NO_OPTION,
                                                             JOptionPane.WARNING_MESSAGE);
            if (status != JOptionPane.YES_OPTION) {
                return false;
            }
        }
        return true;
    }

    private byte[] getValidData(BufferedImage image) {
        if (image != null &&
            image.getColorModel() instanceof IndexColorModel &&
            image.getData().getDataBuffer() instanceof DataBufferByte) {
            return ((DataBufferByte) image.getData().getDataBuffer()).getData();
        }
        return null;
    }

    protected int getNumNonEmptyBins() {
        final byte[] data = getValidData(plot.getImage());
        int n = 0;
        if (data != null) {
            int b;
            for (byte aData : data) {
                b = aData & 0xff;
                if (b != 0) {
                    n++;
                }
            }
        }
        return n;
    }

    @Override
    protected String getDataAsText() {
        final BufferedImage image = plot.getImage();
        final Rectangle2D bounds = plot.getImageDataBounds();

        final byte[] data = getValidData(image);
        if (data == null) {
            return null;
        }

        final StringBuffer sb = new StringBuffer(64000);
        final int w = image.getWidth();
        final int h = image.getHeight();

        final RasterDataNode rasterX = getRaster(X_VAR);
        final String nameX = rasterX.getName();
        final double sampleMinX = bounds.getMinX();
        final double sampleMaxX = bounds.getMaxX();

        final RasterDataNode rasterY = getRaster(Y_VAR);
        final String nameY = rasterY.getName();
        final double sampleMinY = bounds.getMinY();
        final double sampleMaxY = bounds.getMaxY();

        sb.append("Product name:\t").append(rasterX.getProduct().getName()).append("\n");
        sb.append("Dataset X name:\t").append(nameX).append("\n");
        sb.append("Dataset Y name:\t").append(nameY).append("\n");
        sb.append('\n');
        sb.append(nameX).append(" minimum:\t").append(sampleMinX).append("\t").append(rasterX.getUnit()).append(
                "\n");
        sb.append(nameX).append(" maximum:\t").append(sampleMaxX).append("\t").append(rasterX.getUnit()).append(
                "\n");
        sb.append(nameX).append(" bin size:\t").append((sampleMaxX - sampleMinX) / w).append("\t").append(
                rasterX.getUnit()).append("\n");
        sb.append(nameX).append(" #bins:\t").append(w).append("\n");
        sb.append('\n');
        sb.append(nameY).append(" minimum:\t").append(sampleMinY).append("\t").append(rasterY.getUnit()).append(
                "\n");
        sb.append(nameY).append(" maximum:\t").append(sampleMaxY).append("\t").append(rasterY.getUnit()).append(
                "\n");
        sb.append(nameY).append(" bin size:\t").append((sampleMaxY - sampleMinY) / h).append("\t").append(
                rasterY.getUnit()).append("\n");
        sb.append(nameY).append(" #bins:\t").append(h).append("\n");
        sb.append('\n');

        sb.append(nameX);
        sb.append('\t');
        sb.append(nameY);
        sb.append('\t');
        sb.append("Bin counts\t(cropped at 255)");
        sb.append('\n');

        int x, y, z;
        double v1, v2;
        for (int i = 0; i < data.length; i++) {
            z = data[i] & 0xff;
            if (z != 0) {

                x = i % w;
                y = h - i / w - 1;

                v1 = sampleMinX + ((x + 0.5) * (sampleMaxX - sampleMinX)) / w;
                v2 = sampleMinY + ((y + 0.5) * (sampleMaxY - sampleMinY)) / h;

                sb.append(v1);
                sb.append('\t');
                sb.append(v2);
                sb.append('\t');
                sb.append(z);
                sb.append('\n');
            }
        }

        return sb.toString();
    }


    @Override
    public void handleLayerContentChanged() {
        computePanel.updateRoiCheckBoxState();
    }
}

