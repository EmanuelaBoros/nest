package org.esa.beam.visat.toolviews.layermanager.editors;

import com.bc.ceres.glayer.Layer;
import com.bc.ceres.swing.TableLayout;

import org.esa.beam.framework.ui.AppContext;
import org.esa.beam.framework.ui.layer.LayerEditor;
import org.esa.beam.visat.toolviews.layermanager.layersrc.shapefile.FeatureLayer;

import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JSlider;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import java.awt.Font;
import java.util.Hashtable;

/**
 * Editor for placemark layers.
 *
 * @author Ralf Quast
 * @author Marco Zühlke
 * @author Marco Peters
 * @version $Revision: 1.4 $ $Date: 2009-12-22 17:30:01 $
 * @since BEAM 4.6
 */
public class FeatureLayerEditor implements LayerEditor {

    private FeatureLayer currentLayer;
    private JSlider polyFillTransparency;
    private JSlider polyStrokeTransparency;
    private JSlider textTransparency;


    @Override
    public JComponent createControl(AppContext appContext, Layer layer) {
        currentLayer = (FeatureLayer) layer;
        Hashtable<Integer, JLabel> sliderLabelTable = new Hashtable<Integer, JLabel>();
        sliderLabelTable.put(0, createSliderLabel("0%"));
        sliderLabelTable.put(127, createSliderLabel("50%"));
        sliderLabelTable.put(255, createSliderLabel("100%"));

        TableLayout tableLayout = new TableLayout(2);
        tableLayout.setTableFill(TableLayout.Fill.HORIZONTAL);
        tableLayout.setTableAnchor(TableLayout.Anchor.NORTHWEST);
        tableLayout.setColumnWeightX(0, 0.4);
        tableLayout.setColumnWeightX(1, 0.6);
        tableLayout.setRowWeightY(3, 1.0);
        tableLayout.setTablePadding(4, 4);
        JPanel control = new JPanel(tableLayout);

        JLabel fillLabel = new JLabel("Fill transparency:");
        control.add(fillLabel);
        polyFillTransparency = new JSlider(0, 255, 255);
        polyFillTransparency.setToolTipText("Set the opacity of fillings");
        polyFillTransparency.setLabelTable(sliderLabelTable);
        polyFillTransparency.setPaintLabels(true);
        polyFillTransparency.addChangeListener(new ChangeListener() {
            @Override
            public void stateChanged(ChangeEvent e) {
                currentLayer.setPolyFillOpacity(1.0 - polyFillTransparency.getValue() / 255.0);

            }
        });
        control.add(polyFillTransparency);

        JLabel lineLabel = new JLabel("Line transparency:");
        control.add(lineLabel);
        polyStrokeTransparency = new JSlider(0, 255, 255);
        polyStrokeTransparency.setToolTipText("Set the transparency of lines");
        polyStrokeTransparency.setLabelTable(sliderLabelTable);
        polyStrokeTransparency.setPaintLabels(true);
        polyStrokeTransparency.addChangeListener(new ChangeListener() {
            @Override
            public void stateChanged(ChangeEvent e) {
                currentLayer.setPolyStrokeOpacity(1.0 - polyStrokeTransparency.getValue() / 255.0);

            }
        });
        control.add(polyStrokeTransparency);

        JLabel labelLabel = new JLabel("Label transparency:");
        control.add(labelLabel);
        textTransparency = new JSlider(0, 255, 255);
        textTransparency.setToolTipText("Set the transparency of labels");
        textTransparency.setLabelTable(sliderLabelTable);
        textTransparency.setPaintLabels(true);
        textTransparency.addChangeListener(new ChangeListener() {
            @Override
            public void stateChanged(ChangeEvent e) {
                currentLayer.setTextOpacity(1.0 - textTransparency.getValue() / 255.0);

            }
        });
        control.add(textTransparency);
        control.add(new JPanel()); // filler
        return control;
    }

    private JLabel createSliderLabel(String text) {
        JLabel label = new JLabel(text);
        Font oldFont = label.getFont();
        Font newFont = oldFont.deriveFont(oldFont.getSize2D() * 0.85f);
        label.setFont(newFont);
        return label;
    }

    @Override
    public void updateControl() {
        polyFillTransparency.setValue((int) Math.round((1.0 - currentLayer.getPolyFillOpacity()) * 255));
        polyStrokeTransparency.setValue((int) Math.round((1.0 - currentLayer.getPolyStrokeOpacity()) * 255));
        textTransparency.setValue((int) Math.round((1.0 - currentLayer.getTextOpacity()) * 255));
    }
}