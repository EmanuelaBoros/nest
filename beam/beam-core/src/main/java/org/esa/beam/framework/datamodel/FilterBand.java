/*
 * $Id: FilterBand.java,v 1.1 2009-04-28 14:39:32 lveci Exp $
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

import com.bc.ceres.core.ProgressMonitor;
import org.esa.beam.util.Guardian;

import java.awt.Rectangle;
import java.awt.image.Raster;
import java.awt.image.RenderedImage;
import java.io.IOException;

/**
 * Represents a band that generates its data by using another band as input and performs some kind of operation on this input.
 * <p/>
 * <p><i>Note that this class is not yet public API and may change in future releases.</i></p>
 *
 * @author Norman Fomferra
 * @version $Revision: 1.1 $ $Date: 2009-04-28 14:39:32 $
 */
public abstract class FilterBand extends Band {

    private RasterDataNode _source;

    protected FilterBand(String name, int dataType, int width, int height, RasterDataNode source) {
        super(name, dataType, width, height);
        Guardian.assertNotNull("source", source);
        _source = source;
        setSynthetic(true);
    }

    public RasterDataNode getSource() {
        return _source;
    }

    @Override
    public void readRasterData(int offsetX, int offsetY, int width, int height, ProductData rasterData,
                               ProgressMonitor pm) throws IOException {
        pm.beginTask("Reading data...", 2);
        try {
            final RenderedImage sourceImage = getSourceImage();
            final Raster data = sourceImage.getData(new Rectangle(offsetX, offsetY, width, height));
            pm.worked(1);
            data.getDataElements(offsetX, offsetY, width, height, rasterData.getElems());
            pm.worked(1);
        } finally {
            pm.done();
        }
    }

    @Override
    public void dispose() {
        _source = null;
        super.dispose();
    }
}
