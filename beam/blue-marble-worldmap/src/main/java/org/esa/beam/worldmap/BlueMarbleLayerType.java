package org.esa.beam.worldmap;

import com.bc.ceres.binding.PropertyContainer;
import com.bc.ceres.binding.Property;
import com.bc.ceres.glayer.Layer;
import com.bc.ceres.glayer.LayerContext;
import com.bc.ceres.glayer.support.ImageLayer;
import com.bc.ceres.glevel.MultiLevelSource;
import org.esa.beam.glevel.TiledFileMultiLevelSource;
import org.geotools.referencing.AbstractIdentifiedObject;
import org.geotools.referencing.crs.DefaultGeographicCRS;

import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.net.URL;

/**
 * @author Marco Peters
 * @version $Revision: 1.4 $ $Date: 2009-11-04 17:04:32 $
 * @since BEAM 4.6
 */
public class BlueMarbleLayerType extends ImageLayer.Type {

    private static final String WORLD_IMAGE_DIR_PROPERTY_NAME = "org.esa.beam.worldImageDir";
    private volatile MultiLevelSource multiLevelSource;
    private static final String WORLD_MAP_LAYER_NAME = "World Map (NASA Blue Marble)";

    @Override
    public String getName() {
        return "NASA Blue Marble";
    }

    @Override
    public boolean isValidFor(LayerContext ctx) {
        if (ctx.getCoordinateReferenceSystem() instanceof AbstractIdentifiedObject) {
            AbstractIdentifiedObject crs = (AbstractIdentifiedObject) ctx.getCoordinateReferenceSystem();
            return DefaultGeographicCRS.WGS84.equals(crs, false);
        }
        return false;
    }

    @Override
    public Layer createLayer(LayerContext ctx, PropertyContainer configuration) {
        if (multiLevelSource == null) {
            synchronized (this) {
                if (multiLevelSource == null) {
                    multiLevelSource = createMultiLevelSource();
                }
            }
        }
        for (final Property model : super.createLayerConfig(ctx).getProperties()) {
            if (configuration.getProperty(model.getDescriptor().getName()) == null) {
                configuration.addProperty(model);
            }
        }
        configuration.setValue(ImageLayer.PROPERTY_NAME_MULTI_LEVEL_SOURCE, multiLevelSource);

        final ImageLayer layer = new ImageLayer(this, multiLevelSource, configuration);
        layer.setName(WORLD_MAP_LAYER_NAME);
        layer.setVisible(true);

        return layer;
    }

    @Override
    public PropertyContainer createLayerConfig(LayerContext ctx) {
        return new PropertyContainer();
    }

    private static MultiLevelSource createMultiLevelSource() {
        String dirPath = System.getProperty(WORLD_IMAGE_DIR_PROPERTY_NAME);
        if (dirPath == null || dirPath.isEmpty()) {
            dirPath = getDirPathFromModule();
        }
        if (dirPath == null) {
            throw new IllegalStateException("World image directory not found.");
        }
        final MultiLevelSource multiLevelSource;
        try {
            multiLevelSource = TiledFileMultiLevelSource.create(new File(dirPath), false);
        } catch (IOException e) {
            throw new IllegalStateException(e);
        }
        return multiLevelSource;
    }

    private static String getDirPathFromModule() {
        final URL resource = BlueMarbleLayerType.class.getResource("image.properties");
        try {
            return new File(resource.toURI()).getParent();
        } catch (URISyntaxException e) {
            e.printStackTrace();
        }
        return null;
    }

}
