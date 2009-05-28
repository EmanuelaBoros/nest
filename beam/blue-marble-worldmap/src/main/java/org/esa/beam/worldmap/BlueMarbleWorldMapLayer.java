package org.esa.beam.worldmap;

import com.bc.ceres.binding.ValueContainer;
import com.bc.ceres.glayer.LayerType;
import com.bc.ceres.glayer.support.ImageLayer;

/**
 * Provides a world map layer for the SMOS-Box.
 *
 * @author Marco Peters
 * @version $Revision: 1.2 $ $Date: 2009-05-28 14:17:58 $
 * @since BEAM 4.6
 */
public class BlueMarbleWorldMapLayer extends ImageLayer {

    private static final String WORLD_MAP_LAYER_NAME = "World Map (NASA Blue Marble)";

    BlueMarbleWorldMapLayer(ValueContainer multiLevelSource) {
        super((Type) LayerType.getLayerType(BlueMarbleLayerType.class.getName()), multiLevelSource);

        setName(WORLD_MAP_LAYER_NAME);
        setVisible(true);
    }
}
