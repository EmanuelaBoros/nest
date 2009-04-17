package com.bc.ceres.glayer;

import static org.junit.Assert.*;
import org.junit.Test;
import com.bc.ceres.binding.ValidationException;

public class LayerTest {

    @Test
    public void testType() {
        Layer layer1 = new CollectionLayer();
        Layer layer2 = new CollectionLayer();

        assertNotNull(layer1.getLayerType());
        assertSame(layer1.getLayerType(), layer2.getLayerType());
        assertSame(layer1.getLayerType(), LayerType.getLayerType(CollectionLayer.Type.class.getName()));
    }

    @Test
    public void testDefaults() throws ValidationException {
        Layer layer = new CollectionLayer();
        LayerType layerType = layer.getLayerType();
        assertNotNull(layerType);
        assertTrue(layerType.createLayer(null, null) != null);
        assertTrue(layerType.isValidFor(null));
        assertNotNull(layer.getId());
        assertEquals("Collection layer", layer.getName());
        assertEquals(true, layer.isVisible());
        assertEquals(1.0, layer.getStyle().getOpacity(), 1.0e-10);
        assertEquals(Composite.SRC_OVER, layer.getStyle().getComposite());

        assertNull(layer.getModelBounds());
    }

    @Test
    public void testPropertyAccess() {
        Layer layer = new CollectionLayer();

        assertEquals(true, layer.isVisible());
        layer.setVisible(false);
        assertEquals(false, layer.isVisible());
        layer.setVisible(true);
        assertEquals(true, layer.isVisible());

        assertEquals("Collection layer", layer.getName());
        layer.setName("Grid");
        assertEquals("Grid", layer.getName());
        layer.setName("Earth grid");
        assertEquals("Earth grid", layer.getName());
        try {
            layer.setName(null);
            fail("NullPointerException expected");
        } catch (NullPointerException ignored) {
            // expected
        }

        assertEquals(1.0, layer.getStyle().getOpacity(), 1.0e-10);
        layer.getStyle().setOpacity(0.1);
        assertEquals(0.1, layer.getStyle().getOpacity(), 1.0e-10);
        layer.getStyle().setOpacity(1.0);
        assertEquals(1.0, layer.getStyle().getOpacity(), 1.0e-10);

        assertEquals(Composite.SRC_OVER, layer.getStyle().getComposite());
        layer.getStyle().setComposite(Composite.DST_OUT);
        assertEquals(Composite.DST_OUT, layer.getStyle().getComposite());
    }

    @Test
    public void testPropertyChangeNotification() {
        final Layer layer = new CollectionLayer();
        final TracingLayerListener ll = new TracingLayerListener();
        layer.addListener(ll);

        layer.setName("Grid");
        layer.setVisible(false);
        layer.getStyle().setOpacity(0.5);
        assertEquals("name;visible;opacity;", ll.trace);

        layer.setName("Grid");
        layer.setVisible(false);
        layer.getStyle().setOpacity(0.0);
        assertEquals("name;visible;opacity;opacity;", ll.trace);

        layer.getStyle().setOpacity(0.0);
        layer.setVisible(true);
        layer.setName("Raster");
        assertEquals("name;visible;opacity;opacity;visible;name;", ll.trace);

        layer.getStyle().setComposite(Composite.DST_IN);
        assertEquals("name;visible;opacity;opacity;visible;name;composite;", ll.trace);

        ll.trace = "";
        layer.removeListener(ll);

        layer.getStyle().setOpacity(0.25);
        layer.setVisible(false);
        layer.setName("Graticule");
        layer.getStyle().setComposite(Composite.SRC_OUT);
        assertEquals("", ll.trace);
    }

}