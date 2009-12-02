/*
 * $Id: MaskCollectionLayer.java,v 1.2 2009-12-02 16:52:11 lveci Exp $
 *
 * Copyright (C) 2008 by Brockmann Consult (info@brockmann-consult.de)
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation. This program is distributed in the hope it will
 * be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */
package org.esa.beam.glayer;

import com.bc.ceres.binding.PropertyContainer;
import com.bc.ceres.core.Assert;
import com.bc.ceres.glayer.CollectionLayer;
import com.bc.ceres.glayer.Layer;
import com.bc.ceres.glayer.support.AbstractLayerListener;
import com.bc.ceres.glayer.support.ImageLayer;
import org.esa.beam.framework.datamodel.Mask;
import org.esa.beam.framework.datamodel.Product;
import org.esa.beam.framework.datamodel.ProductNode;
import org.esa.beam.framework.datamodel.ProductNodeEvent;
import org.esa.beam.framework.datamodel.ProductNodeGroup;
import org.esa.beam.framework.datamodel.ProductNodeListener;
import org.esa.beam.framework.datamodel.RasterDataNode;

import java.beans.PropertyChangeEvent;
import java.util.HashMap;
import java.util.HashSet;


public class MaskCollectionLayer extends CollectionLayer {

    public static final String ID = MaskCollectionLayer.class.getName();

    private final ProductNodeListener maskPNL;
    private RasterDataNode raster;

    public MaskCollectionLayer(MaskCollectionLayerType layerType,
                               RasterDataNode raster,
                               PropertyContainer configuration) {
        super(layerType, configuration, "Masks");
        Assert.notNull(raster, "raster");
        this.raster = raster;
        this.maskPNL = new MaskPNL();
        setId(ID);
        getProduct().addProductNodeListener(maskPNL);
        addListener(new VisibilityLL());
    }

    @Override
    public void disposeLayer() {
        if (raster != null) {
            getProduct().removeProductNodeListener(maskPNL);
            raster = null;
        }
    }

    private Product getProduct() {
        return raster.getProduct();
    }

    private RasterDataNode getRaster() {
        return raster;
    }

    private Layer createLayer(final Mask mask) {
        return MaskLayerType.createLayer(getRaster(), mask);
    }

    private ImageLayer getMaskLayer(Mask mask) {
        for (Layer layer : getChildren()) {
            if (layer instanceof ImageLayer) {
                final Object value = layer.getConfiguration().getValue(MaskLayerType.PROPERTY_NAME_MASK);
                if (mask == value) {
                    return (ImageLayer) layer;
                }
            }
        }
        return null;
    }

    synchronized void updateChildren() {

        // Collect all current mask layers
        HashMap<Mask, Layer> currentLayers = new HashMap<Mask, Layer>();
        for (Layer maskLayer : getChildren()) {
            Mask mask = (Mask) maskLayer.getConfiguration().getValue("mask");
            currentLayers.put(mask, maskLayer);
        }

        // Allign mask layers with available masks
        Mask[] availableMasks = raster.getProduct().getMaskGroup().toArray(new Mask[0]);
        HashSet<Layer> unusedLayers = new HashSet<Layer>(getChildren());
        for (Mask availableMask : availableMasks) {
            Layer layer = currentLayers.get(availableMask);
            if (layer != null) {
                unusedLayers.remove(layer);
            } else {
                layer = createLayer(availableMask);
                getChildren().add(layer);
            }
            layer.setVisible(raster.getOverlayMaskGroup().contains(availableMask));
        }

        // Remove unused layers
        for (Layer layer : unusedLayers) {
            layer.dispose();
            getChildren().remove(layer);
        }
    }

    public class MaskPNL implements ProductNodeListener {

        @Override
        public synchronized void nodeChanged(ProductNodeEvent event) {
            final ProductNode sourceNode = event.getSourceNode();
            if (sourceNode instanceof Mask) {
                final Mask mask = (Mask) sourceNode;
                final ImageLayer maskLayer = getMaskLayer(mask);
                if (maskLayer != null) {
                    maskLayer.regenerate();
                }
            }
        }

        @Override
        public void nodeDataChanged(ProductNodeEvent event) {
            nodeChanged(event);
        }

        @Override
        public void nodeAdded(ProductNodeEvent event) {
            updateChildren();
        }

        @Override
        public void nodeRemoved(ProductNodeEvent event) {
            updateChildren();
        }
    }

    private class VisibilityLL extends AbstractLayerListener {
        @Override
        public void handleLayerPropertyChanged(Layer layer, PropertyChangeEvent event) {
            if ("visible".equals(event.getPropertyName())) {
                final Object value = layer.getConfiguration().getValue("mask");
                if (value instanceof Mask) {
                    Mask mask = (Mask) value;
                    final ProductNodeGroup<Mask> overlayMaskGroup = getRaster().getOverlayMaskGroup();
                    if (layer.isVisible()) {
                        if (!overlayMaskGroup.contains(mask)) {
                            overlayMaskGroup.add(mask);
                        }
                    } else {
                        overlayMaskGroup.remove(mask);
                    }
                }
            }
        }
    }
}