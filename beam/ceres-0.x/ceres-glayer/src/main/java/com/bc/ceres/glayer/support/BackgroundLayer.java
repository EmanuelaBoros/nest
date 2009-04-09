package com.bc.ceres.glayer.support;

import com.bc.ceres.glayer.Layer;
import com.bc.ceres.glayer.LayerContext;
import com.bc.ceres.glayer.LayerType;
import com.bc.ceres.grender.Rendering;

import java.awt.Graphics2D;
import java.awt.Paint;
import java.awt.Rectangle;
import java.util.HashMap;
import java.util.Map;

/**
 * A background layer is used to draw a background using a unique {@link java.awt.Paint}.
 *
 * @author Norman Fomferra
 */
public class BackgroundLayer extends Layer {
    private static final Type LAYER_TYPE = (Type) LayerType.getLayerType(Type.class.getName());
    
    public BackgroundLayer(Paint paint) {
        this(LAYER_TYPE, paint);
    }
    
    protected BackgroundLayer(Type type, Paint paint) {
        super(type);
        getStyle().setProperty("paint", paint);
    }

    public Paint getPaint() {
        return (Paint) getStyle().getProperty("paint");
    }

    public void setPaint(Paint paint) {
        getStyle().setProperty("paint", paint);
    }

    @Override
    protected void renderLayer(Rendering rendering) {
        final Graphics2D g = rendering.getGraphics();
        Paint oldPaint = g.getPaint();
        g.setPaint(getPaint());
        Rectangle bounds = g.getClipBounds();
        g.fillRect(bounds.x, bounds.y, bounds.width, bounds.height);
        g.setPaint(oldPaint);
    }
    
    public static class Type extends LayerType {
        
        @Override
        public String getName() {
            return "Background Layer";
        }

        @Override
        public boolean isValidFor(LayerContext ctx) {
            return true;
        }

        @Override
        public Map<String, Object> createConfiguration(LayerContext ctx, Layer layer) {
            final HashMap<String, Object> configuration = new HashMap<String, Object>();
            configuration.put("paint", ((BackgroundLayer) layer).getPaint());
            return configuration;
        }

        @Override
        public Layer createLayer(LayerContext ctx, Map<String, Object> configuration) {
            Paint paint = (Paint) configuration.get("paint");
            return new BackgroundLayer(paint);
        }
    }
}