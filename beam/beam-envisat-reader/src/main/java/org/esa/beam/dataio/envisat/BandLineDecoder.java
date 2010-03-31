/*
 * $Id: BandLineDecoder.java,v 1.2 2010-03-31 13:59:56 lveci Exp $
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
package org.esa.beam.dataio.envisat;

/**
 * A decoder which knows how to decode the records of measurement datasets (MDS) in order to create the samples of the
 * related geo-physical (spectral) band raster.
 *
 * @author Norman Fomferra
 * @version $Revision: 1.2 $ $Date: 2010-03-31 13:59:56 $
 */
public interface BandLineDecoder {


    /**
     * Computes the samples of a raster line.
     *
     * @param sourceArray the array containing the values for each x
     * @param sourceMinX  the minimum index for x
     * @param sourceMaxX  the maximum index for x
     * @param sourceStepX the period or step for x
     * @param rasterArray the destination raster array
     * @param rasterPos   the absolute position (offset) within the destination raster array
     * @param rasterIncr  the increment to be used (<code>-1</code> or <code>+1</code>)
     */
    void computeLine(Object sourceArray,
                     int sourceMinX,
                     int sourceMaxX,
                     int sourceStepX,
                     Object rasterArray,
                     int rasterPos,
                     int rasterIncr);
}
