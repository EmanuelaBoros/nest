/*
 * $Id: BandImagesTest.java,v 1.1 2009-11-04 17:04:32 lveci Exp $
 *
 * Copyright (C) 2002 by Brockmann Consult (info@brockmann-consult.de)
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
package org.esa.beam.framework.datamodel;

import org.esa.beam.jai.ImageManager;
import static org.junit.Assert.*;
import org.junit.Test;

import javax.media.jai.PlanarImage;
import javax.media.jai.TiledImage;
import java.awt.image.ColorModel;
import java.awt.image.DataBuffer;
import java.awt.image.PixelInterleavedSampleModel;
import java.awt.image.SampleModel;

public class BandImagesTest {

    @Test
    public void testGeophysicalSameAsSourceImage() {
        Product p = new Product("p", "pt", 256, 128);
        Band band;

        band = addBand(p, "b1", ProductData.TYPE_UINT16, null, 1.0);
        assertNotNull(band.getSourceImage());
        assertEquals(DataBuffer.TYPE_USHORT, band.getSourceImage().getSampleModel().getDataType());
        assertSame(band.getGeophysicalImage(), band.getSourceImage());
        assertNull(band.getValidMaskImage());

        band = addBand(p, "b2", ProductData.TYPE_FLOAT32, null, 1.0);
        assertNotNull(band.getSourceImage());
        assertEquals(DataBuffer.TYPE_FLOAT, band.getSourceImage().getSampleModel().getDataType());
        assertSame(band.getGeophysicalImage(), band.getSourceImage());
        assertNull(band.getValidMaskImage());
    }

    @Test
    public void testGeophysicalDifferentFromSourceImage() {
        Product p = new Product("p", "pt", 256, 128);
        Band band;

        band = addBand(p, "b1", ProductData.TYPE_UINT16, null, 1.1);
        assertNotNull(band.getSourceImage());
        assertEquals(DataBuffer.TYPE_USHORT, band.getSourceImage().getSampleModel().getDataType());
        assertNotSame(band.getGeophysicalImage(), band.getSourceImage());
        assertEquals(DataBuffer.TYPE_FLOAT, band.getGeophysicalImage().getSampleModel().getDataType());
        assertNotNull(band.getGeophysicalImage());
        assertNull(band.getValidMaskImage());

        band = addBand(p, "b2", ProductData.TYPE_FLOAT32, null, 1.1);
        assertNotNull(band.getSourceImage());
        assertEquals(DataBuffer.TYPE_FLOAT, band.getSourceImage().getSampleModel().getDataType());
        assertNotSame(band.getGeophysicalImage(), band.getSourceImage());
        assertEquals(DataBuffer.TYPE_DOUBLE, band.getGeophysicalImage().getSampleModel().getDataType());
        assertNotNull(band.getGeophysicalImage());
        assertNull(band.getValidMaskImage());
    }

    @Test
    public void testValidMaskImage() {
        Product p = new Product("p", "pt", 256, 128);
        Band band;
        band = addBand(p, "b1", ProductData.TYPE_UINT16, -999.0, 1.0);
        assertNotNull(band.getSourceImage());
        assertEquals(DataBuffer.TYPE_USHORT, band.getSourceImage().getSampleModel().getDataType());
        assertSame(band.getGeophysicalImage(), band.getSourceImage());
        assertNotNull(band.getValidMaskImage());
        assertEquals(DataBuffer.TYPE_BYTE, band.getValidMaskImage().getSampleModel().getDataType());
    }

    private static Band addBand(Product p, String name, int dataType, Double noDataValue, double scalingFactor) {
        final int width = p.getSceneRasterWidth();
        final int height = p.getSceneRasterHeight();
        final SampleModel sm = new PixelInterleavedSampleModel(ImageManager.getDataBufferType(dataType), width, height, 1, width, new int[]{0});
        final ColorModel cm = PlanarImage.createColorModel(sm);
        final TiledImage sourceImage = new TiledImage(0, 0, width, height, 0, 0, sm, cm);


        final Band band = new Band(name, dataType, width, height);
        band.setScalingFactor(scalingFactor);
        if (noDataValue != null) {
            band.setNoDataValueUsed(true);
            band.setNoDataValue(noDataValue);
        }
        band.setSourceImage(sourceImage);

        p.addBand(band);

        return band;
    }
}