/*
 * $Id: RoiLayerTypeTest.java,v 1.3 2010-02-12 14:42:16 lveci Exp $
 *
 * Copyright (C) 2009 by Brockmann Consult (info@brockmann-consult.de)
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

import com.bc.ceres.binding.PropertySet;
import com.bc.ceres.glayer.Layer;
import com.bc.ceres.glayer.LayerTypeRegistry;
import com.bc.ceres.glayer.support.ImageLayer;
import static junit.framework.Assert.assertSame;
import org.esa.beam.framework.datamodel.Band;
import org.esa.beam.framework.datamodel.Mask;
import org.esa.beam.framework.datamodel.Product;
import org.esa.beam.framework.datamodel.ProductData;
import org.esa.beam.framework.datamodel.RasterDataNode;
import org.esa.beam.framework.datamodel.VirtualBand;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;
import org.junit.Test;

import java.awt.Color;


public class RoiLayerTypeTest extends LayerTypeTest {

    public RoiLayerTypeTest() {
        super(RoiLayerType.class);
    }

    @Test
    public void testConfigurationTemplate() {
        final PropertySet template = getLayerType().createLayerConfig(null);

        assertNotNull(template);
        ensurePropertyIsDeclaredButNotDefined(template, RoiLayerType.PROPERTY_NAME_RASTER, RasterDataNode.class);

        ensurePropertyIsDefined(template, RoiLayerType.PROPERTY_NAME_COLOR, Color.class);
        ensurePropertyIsDefined(template, RoiLayerType.PROPERTY_NAME_TRANSPARENCY, Double.class);
    }
    
    @Test
    public void testCreateLayer() {
        final Product product = new Product("N", "T", 10, 10);
        final Band raster = new VirtualBand("A", ProductData.TYPE_INT32, 10, 10, "42");
        product.addBand(raster);
        final Mask mask = new Mask("A_roi", 10, 10, Mask.BandMathType.INSTANCE);
        Mask.BandMathType.setExpression(mask, "A == 42");
        product.getMaskGroup().add(mask);

        final PropertySet config = getLayerType().createLayerConfig(null);
        config.setValue(RoiLayerType.PROPERTY_NAME_RASTER, raster);

        final Layer layer = getLayerType().createLayer(null, config);
        assertNotNull(layer);
        
        MaskLayerType maskLayerType = LayerTypeRegistry.getLayerType(MaskLayerType.class);
        assertSame(maskLayerType, layer.getLayerType());
        assertTrue(layer instanceof ImageLayer);
        ImageLayer imageLayer = (ImageLayer) layer;
        assertNotNull(imageLayer.getMultiLevelSource());
    }
}
