package org.esa.beam.visat.toolviews.layermanager;

import com.bc.ceres.glayer.Layer;

import org.esa.beam.framework.ui.AppContext;

import javax.swing.JComponent;

/**
 * An editor for a specific layer type.
 *
 * @author Norman Fomferra
 * @version $Revision: 1.1 $ $Date: 2009-04-27 13:08:25 $
 * @since BEAM 4.6
 */
public interface LayerEditor {

    /**
     * Creates the control for the user interface which is displayed
     * in the Layer Editor Toolview.
     * 
     * @param appContext the application context 
     * @param layer The layer to create the control for.
     * @return The control.
     */
    JComponent createControl(AppContext appContext, Layer layer);

    /**
     * It is called whenever the control must be updated.
     */
    void updateControl();
}
