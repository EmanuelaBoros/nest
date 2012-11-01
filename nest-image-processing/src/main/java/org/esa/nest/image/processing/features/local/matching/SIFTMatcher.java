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
import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.esa.beam.dataio.envisat.EnvisatProductReaderPlugIn;
import org.esa.beam.framework.datamodel.Band;
import org.esa.beam.framework.gpf.OperatorSpi;
import org.esa.beam.framework.gpf.Tile;
import org.esa.nest.image.processing.features.local.matching.SARTemplateMatcher.Mode;
import org.openimaj.feature.local.list.LocalFeatureList;
import org.openimaj.feature.local.matcher.FastBasicKeypointMatcher;
import org.openimaj.feature.local.matcher.MatchingUtilities;
import org.openimaj.feature.local.matcher.consistent.ConsistentLocalFeatureMatcher2d;
import org.openimaj.image.DisplayUtilities;
import org.openimaj.image.FImage;
import org.openimaj.image.ImageUtilities;
import org.openimaj.image.MBFImage;
import org.openimaj.image.colour.RGBColour;
import org.openimaj.image.feature.local.engine.DoGSIFTEngine;
import org.openimaj.image.feature.local.keypoints.Keypoint;
import org.openimaj.math.geometry.point.Point2d;
import org.openimaj.math.geometry.transforms.HomographyModel;
import org.openimaj.math.model.fit.RANSAC;
import org.openimaj.util.pair.Pair;

/**
 *
 * @author Emanuela
 */
public class SIFTMatcher implements TemplateMatcherEnforcement {

    DoGSIFTEngine engine;

    @Override
    public float computeMatchScore(final Band sourceBand, final Tile sourceRaster,
            final Band targetBand, final Tile targetTile, Mode mode) {

        engine = new DoGSIFTEngine();

        RenderedImage fullRenderedImage = sourceBand.getSourceImage().getImage(0);
        BufferedImage fullBufferedImage =
                new BufferedImage(sourceBand.getSceneRasterWidth(),
                sourceBand.getSceneRasterHeight(),
                BufferedImage.TYPE_USHORT_GRAY);
        fullBufferedImage.setData(fullRenderedImage.getData());

        final Rectangle sourceRectangle = sourceRaster.getRectangle();
        BufferedImage imageBuffer = fullBufferedImage.getSubimage(
                sourceRectangle.x, sourceRectangle.y,
                sourceRectangle.width, sourceRectangle.height);
        FImage image = ImageUtilities.createFImage(imageBuffer);

        Rectangle templateRectangle = targetTile.getRectangle();
        BufferedImage templateBuffer = fullBufferedImage.getSubimage(
                templateRectangle.x, templateRectangle.y,
                templateRectangle.width, templateRectangle.height);

        FImage template = ImageUtilities.createFImage(templateBuffer);

        LocalFeatureList<Keypoint> k1 = engine.findFeatures(image);
        LocalFeatureList<Keypoint> k2 = engine.findFeatures(template);

        ConsistentLocalFeatureMatcher2d<Keypoint> matcher =
                new ConsistentLocalFeatureMatcher2d<Keypoint>(
                new FastBasicKeypointMatcher<Keypoint>(8));
        HomographyModel model = new HomographyModel(8);
        RANSAC<Point2d, Point2d> modelFitting = new RANSAC<Point2d, Point2d>(
                model, 1600, new RANSAC.BestFitStoppingCondition(), true);
        matcher.setFittingModel(modelFitting);

        matcher.setModelFeatures(k1);
        matcher.findMatches(k2);

        List<Pair<Keypoint>> matches = matcher.getMatches();

        MBFImage inp1MBF = image.toRGB();
        MBFImage inp2MBF = template.toRGB();

        try {
            ImageUtilities.write(MatchingUtilities.drawMatches(inp1MBF, inp2MBF, matches, RGBColour.RED),
                    new File("D:/cat-sift.png"));
        } catch (IOException ex) {
            Logger.getLogger(SIFTMatcher.class.getName()).log(Level.SEVERE, null, ex);
        }

        float fScore = 0f;
        for (Pair<Keypoint> p : matcher.getMatches()) {
            double accum = 0d;
            byte[] v1 = p.firstObject().ivec;
            byte[] v2 = p.secondObject().ivec;
            for (int i = 0; i < v1.length; i++) {
                double v1i = ((double) v1[i]);
                double v2i = ((double) v2[i]);
                accum += (v1i - v2i) * (v1i - v2i);
            }
            fScore += Math.sqrt(accum);
        }

        return fScore;
    }
}
