package org.esa.beam.visat.toolviews.stat;

import com.bc.ceres.core.ProgressMonitor;
import com.bc.ceres.swing.TableLayout;
import com.bc.ceres.swing.progress.ProgressMonitorSwingWorker;
import org.esa.beam.framework.datamodel.RasterDataNode;
import org.esa.beam.framework.datamodel.Stx;
import org.esa.beam.framework.param.ParamChangeEvent;
import org.esa.beam.framework.param.ParamChangeListener;
import org.esa.beam.framework.param.ParamGroup;
import org.esa.beam.framework.param.Parameter;
import org.esa.beam.framework.ui.GridBagUtils;
import org.esa.beam.framework.ui.application.ToolView;
import org.esa.beam.util.math.MathUtils;
import org.esa.beam.util.math.Range;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.LogarithmicAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYBarRenderer;
import org.jfree.data.xy.XIntervalSeries;
import org.jfree.data.xy.XIntervalSeriesCollection;
import org.jfree.ui.RectangleInsets;

import javax.swing.BorderFactory;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import java.awt.BorderLayout;
import java.awt.Container;
import java.awt.GridBagConstraints;
import java.awt.image.RenderedImage;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

/**
 * A pane within the statistcs window which displays a histogram.
 *
 * @author Marco Peters
 * @version $Revision: 1.1 $ $Date: 2009-04-27 13:08:26 $
 */
class HistogramPanel extends PagePanel {

    private static final String NO_DATA_MESSAGE = "No histogram computed yet.\n" +
                                                  "TIP: To zoom within the chart draw a rectangle\n" +
                                                  "with the mouse or use the context menu.";  /*I18N*/
    private static final String CHART_TITLE = "Histogram";
    private static final String TITLE_PREFIX = CHART_TITLE;    /*I18N*/


    private Parameter numBinsParam;
    private Parameter autoMinMaxEnabledParam;
    private Parameter histoMinParam;
    private Parameter histoMaxParam;
    private ChartPanel histogramDisplay;
    private boolean histogramComputing;
    private ComputePanel computePanel;
    private XIntervalSeriesCollection dataset;

    private JFreeChart chart;
    private Stx stx;


    public HistogramPanel(final ToolView parentDialog, String helpID) {
        super(parentDialog, helpID);
    }

    @Override
    protected String getTitlePrefix() {
        return TITLE_PREFIX;
    }

    public JFreeChart getChart() {
        return chart;
    }

    @Override
    protected void initContent() {
        initParameters();
        createUI();
        updateContent();
    }

    @Override
    protected void updateContent() {
        if (computePanel != null) {
            computePanel.setRaster(getRaster());
            setRaster(getRaster());
            if (!(Boolean) autoMinMaxEnabledParam.getValue()) {
                return;
            }
            chart.setTitle(getRaster() != null ? CHART_TITLE + " for " + getRaster().getName() : CHART_TITLE);
            if (getRaster() != null && getRaster().isLog10Scaled()) {
                ValueAxis oldDomainAxis = chart.getXYPlot().getDomainAxis();
                if (!(oldDomainAxis instanceof LogarithmicAxis)) {
                    LogarithmicAxis logAxisX = new LogarithmicAxis("Values");
                    logAxisX.setAllowNegativesFlag(true);
                    chart.getXYPlot().setDomainAxis(logAxisX);
                }
            } else {
                ValueAxis oldDomainAxis = chart.getXYPlot().getDomainAxis();
                if (oldDomainAxis instanceof LogarithmicAxis) {
                    NumberAxis xAxis = new NumberAxis("Values");
                    chart.getXYPlot().setDomainAxis(xAxis);
                }
            }
            histoMinParam.setDefaultValue();
            histoMaxParam.setDefaultValue();

            if(getRaster() != null) {
                final Stx stat = getRaster().getStx();
                if(stat != null) {
                    setStx(stat);
                }
            }
        }
    }

    @Override
    protected boolean mustUpdateContent() {
        return isRasterChanged();
    }

    private void initParameters() {
        ParamGroup paramGroup = new ParamGroup();

        numBinsParam = new Parameter("histo.numBins", Stx.DEFAULT_BIN_COUNT);
        numBinsParam.getProperties().setLabel("#Bins:");   /*I18N*/
        numBinsParam.getProperties().setDescription("Set the number of bins in the histogram");    /*I18N*/
        numBinsParam.getProperties().setMinValue(2);
        numBinsParam.getProperties().setMaxValue(2000);
        numBinsParam.getProperties().setNumCols(5);
        paramGroup.addParameter(numBinsParam);

        autoMinMaxEnabledParam = new Parameter("histo.autoMinMax", Boolean.TRUE);
        autoMinMaxEnabledParam.getProperties().setLabel("Auto min/max");   /*I18N*/
        autoMinMaxEnabledParam.getProperties().setDescription("Automatically detect min/max"); /*I18N*/
        paramGroup.addParameter(autoMinMaxEnabledParam);

        histoMinParam = new Parameter("histo.min", 0.0);
        histoMinParam.getProperties().setLabel("Min:");    /*I18N*/
        histoMinParam.getProperties().setDescription("Histogram minimum sample value");    /*I18N*/
        histoMinParam.getProperties().setNumCols(7);
        paramGroup.addParameter(histoMinParam);

        histoMaxParam = new Parameter("histo.max", 100.0);
        histoMaxParam.getProperties().setLabel("Max:");    /*I18N*/
        histoMaxParam.getProperties().setDescription("Histogram maximum sample value");    /*I18N*/
        histoMaxParam.getProperties().setNumCols(7);
        paramGroup.addParameter(histoMaxParam);

        paramGroup.addParamChangeListener(new ParamChangeListener() {

            public void parameterValueChanged(ParamChangeEvent event) {
                updateUIState();
            }
        });
    }

    private void createUI() {
        dataset = new XIntervalSeriesCollection();
        chart = ChartFactory.createHistogram(
                CHART_TITLE,
                "Values",
                "Frequency",
                dataset,
                PlotOrientation.VERTICAL,
                false,  // Legend?
                true,   // tooltips
                false   // url
        );
        final XYPlot xyPlot = chart.getXYPlot();
        XYBarRenderer renderer = (XYBarRenderer) xyPlot.getRenderer();
        renderer.setDrawBarOutline(false);
        renderer.setBaseToolTipGenerator(new XYPlotToolTipGenerator());

        histogramDisplay = createChartPanel(chart);

        final ActionListener actionAll = new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                computeHistogram(false);
            }
        };
        final ActionListener actionROI = new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                computeHistogram(true);
            }
        };
        computePanel = ComputePanel.createComputePane(actionAll, actionROI, getRaster());
        final TableLayout rightPanelLayout = new TableLayout(1);
        final JPanel rightPanel = new JPanel(rightPanelLayout);
        rightPanelLayout.setTableFill(TableLayout.Fill.HORIZONTAL);
        rightPanelLayout.setRowWeightY(3, 1.0);
        rightPanelLayout.setRowAnchor(4, TableLayout.Anchor.EAST);
        rightPanel.add(computePanel);
        rightPanel.add(createOptionsPane());
        rightPanel.add(createChartButtonPanel(histogramDisplay));
        rightPanel.add(new JPanel());   // filler
        final JPanel helpPanel = new JPanel(new BorderLayout());
        helpPanel.add(getHelpButton(), BorderLayout.EAST);
        rightPanel.add(helpPanel);

        add(histogramDisplay, BorderLayout.CENTER);
        add(rightPanel, BorderLayout.EAST);
        updateUIState();
    }

    private ChartPanel createChartPanel(JFreeChart chart) {
        XYPlot plot = chart.getXYPlot();

        plot.setForegroundAlpha(0.85f);
        plot.setNoDataMessage(NO_DATA_MESSAGE);
        plot.setAxisOffset(new RectangleInsets(5, 5, 5, 5));

        ChartPanel chartPanel = new ChartPanel(chart);
        chartPanel.getPopupMenu().add(createCopyDataToClipboardMenuItem());
        return chartPanel;
    }

    private void updateUIState() {
        final double min = ((Number) histoMinParam.getValue()).doubleValue();
        final double max = ((Number) histoMaxParam.getValue()).doubleValue();
        if (!histogramComputing && min > max) {
            histoMinParam.setValue(max, null);
            histoMaxParam.setValue(min, null);
        }
        final boolean autoMinMaxEnabled = (Boolean) autoMinMaxEnabledParam.getValue();
        histoMinParam.setUIEnabled(!autoMinMaxEnabled);
        histoMaxParam.setUIEnabled(!autoMinMaxEnabled);
    }

    private JPanel createOptionsPane() {
        final JPanel optionsPane = GridBagUtils.createPanel();
        final GridBagConstraints gbc = GridBagUtils.createConstraints("anchor=WEST,fill=BOTH");

        GridBagUtils.setAttributes(gbc, "gridwidth=1");

        GridBagUtils.setAttributes(gbc, "gridwidth=1,gridy=0,insets.top=2");
        GridBagUtils.addToPanel(optionsPane, numBinsParam.getEditor().getLabelComponent(), gbc,
                                "gridx=0,weightx=1");
        GridBagUtils.addToPanel(optionsPane, numBinsParam.getEditor().getComponent(), gbc, "gridx=1,weightx=0");

        GridBagUtils.setAttributes(gbc, "gridwidth=2,gridy=1,insets.top=7");
        GridBagUtils.addToPanel(optionsPane, autoMinMaxEnabledParam.getEditor().getComponent(), gbc,
                                "gridwidth=2,gridx=0,weightx=1");

        GridBagUtils.setAttributes(gbc, "gridwidth=1,gridy=2,insets.top=2");
        GridBagUtils.addToPanel(optionsPane, histoMinParam.getEditor().getLabelComponent(), gbc,
                                "gridx=0,weightx=1");
        GridBagUtils.addToPanel(optionsPane, histoMinParam.getEditor().getComponent(), gbc, "gridx=1,weightx=0");

        GridBagUtils.setAttributes(gbc, "gridwidth=1,gridy=3,insets.top=2");
        GridBagUtils.addToPanel(optionsPane, histoMaxParam.getEditor().getLabelComponent(), gbc,
                                "gridx=0,weightx=1");
        GridBagUtils.addToPanel(optionsPane, histoMaxParam.getEditor().getComponent(), gbc, "gridx=1,weightx=0");

        optionsPane.setBorder(BorderFactory.createTitledBorder("X"));

        return optionsPane;
    }

    private void computeHistogram(final boolean useROI) {
        final RenderedImage roiImage;
        if (useROI) {
            roiImage = getRoiImage(getRaster());
        } else {
            roiImage = null;
        }

        final int numBins = ((Number) numBinsParam.getValue()).intValue();
        final boolean autoMinMaxEnabled = getAutoMinMaxEnabled();
        final Range range;
        if (autoMinMaxEnabled) {
            range = null; // auto compute range
        } else {
            final double min = ((Number) histoMinParam.getValue()).doubleValue();
            final double max = ((Number) histoMaxParam.getValue()).doubleValue();
            range = new Range(getRaster().scaleInverse(min), getRaster().scaleInverse(max));
        }

        ProgressMonitorSwingWorker<Stx, Object> swingWorker = new ProgressMonitorSwingWorker<Stx, Object>(this, "Computing Histogram") {
            @Override
            protected Stx doInBackground(ProgressMonitor pm) throws Exception {
                final Stx stx;
                if (roiImage == null && range == null && numBins == Stx.DEFAULT_BIN_COUNT) {
                    stx = getRaster().getStx(true, pm);
                } else if (range == null) {
                    stx = Stx.create(getRaster(), roiImage, numBins, pm);
                } else {
                    stx = Stx.create(getRaster(), roiImage, numBins, range.getMin(), range.getMax(), pm);
                }
                return stx;
            }

            @Override
            public void done() {
                try {
                    Stx stx = get();
                    if (stx.getSampleCount() > 0) {
                        if (autoMinMaxEnabled) {
                            final double min = getRaster().scale(stx.getMin());
                            final double max = getRaster().scale(stx.getMax());
                            final double v = MathUtils.computeRoundFactor(min, max, 4);
                            histogramComputing = true;
                            histoMinParam.setValue(StatisticsUtils.round(min, v), null);
                            histoMaxParam.setValue(StatisticsUtils.round(max, v), null);
                            histogramComputing = false;
                        }
                    } else {
                        JOptionPane.showMessageDialog(getParentComponent(),
                                                      "The ROI is empty or no pixels found between min/max.\n"
                                                              + "A valid histogram could not be computed.",
                                                      CHART_TITLE,
                                                      JOptionPane.WARNING_MESSAGE);
                    }
                    setStx(stx);
                } catch (Exception e) {
                    e.printStackTrace();
                    JOptionPane.showMessageDialog(getParentComponent(),
                                                  "Failed to compute histogram.\nAn internal error occured:\n" + e.getMessage(),
                                                  CHART_TITLE,
                                                  JOptionPane.ERROR_MESSAGE);
                }
            }
        };
        swingWorker.execute();
    }

    private void setStx(Stx stx) {
        this.stx = stx;
        dataset = new XIntervalSeriesCollection();
        if (this.stx != null) {
            final int[] binCounts = this.stx.getHistogramBins();
            final RasterDataNode raster = getRaster();
            final XIntervalSeries series = new XIntervalSeries(raster.getName());
            for (int i = 0; i < binCounts.length; i++) {
                final double xMin = raster.scale(stx.getHistogramBinMin(i));
                final double xMax = raster.scale(stx.getHistogramBinMax(i));
                final double xAvg = (xMin + xMax) * 0.5;
                series.add(xAvg, xMin, xMax, binCounts[i]);
            }
            dataset.addSeries(series);
            chart.getXYPlot().getDomainAxis().setLabel(getAxisLabel(raster));
        }
        chart.getXYPlot().setDataset(dataset);
        chart.fireChartChanged();
    }

    private static String getAxisLabel(RasterDataNode raster) {
        final String unit = raster.getUnit();
        if (unit != null) {
            return raster.getName() + " (" + unit + ")";
        }
        return raster.getName();
    }

    private Container getParentComponent() {
        return getParentDialog().getContext().getPane().getControl();
    }

    private boolean getAutoMinMaxEnabled() {
        return (Boolean) autoMinMaxEnabledParam.getValue();
    }


    @Override
    public String getDataAsText() {
        if (stx == null) {
            return null;
        }

        final int[] binVals = stx.getHistogramBins();
        final int numBins = binVals.length;
        final double min = getRaster().scale(stx.getMin());
        final double max = getRaster().scale(stx.getMax());

        final StringBuilder sb = new StringBuilder(16000);

        sb.append("Product name:\t").append(getRaster().getProduct().getName()).append("\n");
        sb.append("Dataset name:\t").append(getRaster().getName()).append("\n");
        sb.append('\n');
        sb.append("Histogram minimum:\t").append(min).append("\t").append(getRaster().getUnit()).append("\n");
        sb.append("Histogram maximum:\t").append(max).append("\t").append(getRaster().getUnit()).append("\n");
        sb.append("Histogram bin size:\t").append(
                getRaster().isLog10Scaled() ? ("NA\t") : ((max - min) / numBins + "\t") +
                        getRaster().getUnit() + "\n");
        sb.append("Histogram #bins:\t").append(numBins).append("\n");
        sb.append('\n');

        sb.append("Bin center value");
        sb.append('\t');
        sb.append("Bin counts");
        sb.append('\n');

        for (int i = 0; i < numBins; i++) {
            sb.append(min + ((i + 0.5) * (max - min)) / numBins);
            sb.append('\t');
            sb.append(binVals[i]);
            sb.append('\n');
        }

        return sb.toString();
    }

    @Override
    public void handleLayerContentChanged() {
        computePanel.updateRoiCheckBoxState();
    }
}

