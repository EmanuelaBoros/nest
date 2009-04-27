package org.esa.beam.visat.toolviews.layermanager;

import com.bc.ceres.glayer.LayerType;

/**
 * A descriptor for a layer editor.
 *
 * @author Norman Fomferra
 * @version $Revision: 1.1 $ $Date: 2009-04-27 13:08:25 $
 * @since BEAM 4.6
 */
public interface LayerEditorDescriptor {

    /**
     * Gets the {@link LayerType} class for which the {@link LayerEditor} is intended.
     * The corresponding {@code LayerEditor} class can be retrieved by a call to
     * {@link #getLayerEditorClass()}.
     *
     * @return The {@link LayerType} class.
     */
    Class<? extends LayerType> getLayerTypeClass();

    /**
     * Gets the {@link LayerEditor} class, which is intended for the
     * {@link LayerType} class returned by {@link #getLayerTypeClass()}
     *
     * @return The {@link LayerEditor} class
     */
    Class<? extends LayerEditor> getLayerEditorClass();
}
