/*
 * Copyright (C) 2012 by Array Systems Computing Inc. http://www.array.ca
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation; either version 3 of the License, or (at your option)
 * any later version.
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, see http://www.gnu.org/licenses/
 */
package org.esa.nest.image.processing.features.local.matching;

import java.awt.Rectangle;
import java.awt.image.BufferedImage;
import java.awt.image.RenderedImage;
import org.esa.beam.framework.datamodel.Band;
import org.esa.beam.framework.gpf.Tile;
import org.openimaj.image.FImage;
import org.openimaj.image.ImageUtilities;
import org.openimaj.image.analysis.algorithm.TemplateMatcher;

/**
 *
 * @author Emanuela
 */
public class SARTemplateMatcher implements TemplateMatcherEnforcement {

    public enum Mode {

        SUM_SQUARED_DIFFERENCE, NORM_SUM_SQUARED_DIFFERENCE, CORRELATION,
        NORM_CORRELATION, CORRELATION_COEFFICIENT, NORM_CORRELATION_COEFFICIENT
    }

    /**
     * Compute the score at a point as the sum-squared difference between the
     * image and the template with the top-left at the given point. The
     * SARTemplateMatcher will account for the offset to the center of the
     * template internally.
     *
     */
    @Override
    public float computeMatchScore(final Band sourceBand, final Tile sourceRaster,
            final Band targetBand, final Tile targetTile,
            Mode mode) {

        final RenderedImage fullRenderedImage = sourceBand.getSourceImage().getImage(0);
        final BufferedImage fullBufferedImage =
                new BufferedImage(sourceBand.getSceneRasterWidth(),
                sourceBand.getSceneRasterHeight(),
                BufferedImage.TYPE_USHORT_GRAY);
        fullBufferedImage.setData(fullRenderedImage.getData());
        final Rectangle srcTileRectangle = sourceRaster.getRectangle();

        BufferedImage imageBuffer = fullBufferedImage.getSubimage(
                srcTileRectangle.x, srcTileRectangle.y,
                srcTileRectangle.width, srcTileRectangle.height);
        FImage image = ImageUtilities.createFImage(imageBuffer);

//        Rectangle templateRectangle = templateTile.getRectangle();
//
//        BufferedImage templateBuffer = fullBufferedImage.getSubimage(
//                templateRectangle.x, templateRectangle.y,
//                templateRectangle.width, templateRectangle.height);
//        FImage template = ImageUtilities.createFImage(templateBuffer);
//
//        final float[][] imageData = image.pixels;
//        final float[][] templateData = template.pixels;
        float fNorm=0f;
//
//        switch (mode) {
//            case SUM_SQUARED_DIFFERENCE:
//                fNorm = TemplateMatcher.Mode.SUM_SQUARED_DIFFERENCE.
//                        computeMatchScore(imageData, 0, 0,
//                        templateData, 0, 0, srcTileRectangle.width,
//                        srcTileRectangle.height);
//                break;
//            case NORM_SUM_SQUARED_DIFFERENCE:
//                fNorm = TemplateMatcher.Mode.NORM_SUM_SQUARED_DIFFERENCE.
//                        computeMatchScore(imageData, 0, 0,
//                        templateData, 0, 0, srcTileRectangle.width,
//                        srcTileRectangle.height);
//                break;
//            case CORRELATION:
//                fNorm = TemplateMatcher.Mode.CORRELATION.
//                        computeMatchScore(imageData, 0, 0,
//                        templateData, 0, 0, srcTileRectangle.width,
//                        srcTileRectangle.height);
//                break;
//            case NORM_CORRELATION:
//                fNorm = TemplateMatcher.Mode.NORM_CORRELATION.
//                        computeMatchScore(imageData, 0, 0,
//                        templateData, 0, 0, srcTileRectangle.width,
//                        srcTileRectangle.height);
//                break;
//            case CORRELATION_COEFFICIENT:
//                fNorm = TemplateMatcher.Mode.CORRELATION_COEFFICIENT.
//                        computeMatchScore(imageData, 0, 0,
//                        templateData, 0, 0, srcTileRectangle.width,
//                        srcTileRectangle.height);
//                break;
//            case NORM_CORRELATION_COEFFICIENT:
//                fNorm = TemplateMatcher.Mode.CORRELATION_COEFFICIENT.
//                        computeMatchScore(imageData, 0, 0,
//                        templateData, 0, 0, srcTileRectangle.width,
//                        srcTileRectangle.height);
//                break;
//            default:
//                fNorm = 0f;
//        }
        return fNorm;
    }
}