/*
 * $Id: Graticule.java,v 1.1 2009-04-28 14:39:32 lveci Exp $
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

import org.esa.beam.util.Debug;
import org.esa.beam.util.Guardian;
import org.esa.beam.util.ProductUtils;

import java.awt.*;
import java.awt.geom.GeneralPath;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

/**
 * A geometric representation of a geographical grid measured in longitudes and latitudes.
 */
public class Graticule {

    private final GeneralPath[] _linePaths;
    private final TextGlyph[] _textGlyphs;

    private Graticule(GeneralPath[] paths, TextGlyph[] textGlyphs) {
        _linePaths = paths;
        _textGlyphs = textGlyphs;
    }

    public GeneralPath[] getLinePaths() {
        return _linePaths;
    }

    public TextGlyph[] getTextGlyphs() {
        return _textGlyphs;
    }

    /**
     * Creates a graticule for the given product.
     *
     * @param product              the product
     * @param autoDeterminingSteps if true, <code>gridCellSize</code> is used to compute <code>latMajorStep</code>, <code>lonMajorStep</code> for the given product
     * @param gridCellSize         the grid cell size in pixels, ignored if <code>autoDeterminingSteps</code> if false
     * @param latMajorStep         the grid cell size in meridional direction, ignored if <code>autoDeterminingSteps</code> if true
     * @param lonMajorStep         the grid cell size in parallel direction, ignored if <code>autoDeterminingSteps</code> if true
     * @return the graticule or null, if it could not be created
     */
    public static Graticule create(Product product,
                                   boolean autoDeterminingSteps,
                                   int gridCellSize,
                                   float latMajorStep,
                                   float lonMajorStep) {
        Guardian.assertNotNull("product", product);
        final GeoCoding geoCoding = product.getGeoCoding();
        if (geoCoding == null || product.getSceneRasterWidth() < 16 || product.getSceneRasterHeight() < 16) {
            return null;
        }

        if (autoDeterminingSteps) {
            final PixelPos pixelPos1 = new PixelPos(0.5f * product.getSceneRasterWidth(), 0.5f * product.getSceneRasterHeight());
            final PixelPos pixelPos2 = new PixelPos(pixelPos1.x + 1f, pixelPos1.y + 1f);
            final GeoPos geoPos1 = geoCoding.getGeoPos(pixelPos1, null);
            final GeoPos geoPos2 = geoCoding.getGeoPos(pixelPos2, null);
            double deltaLat = Math.abs(geoPos2.lat - geoPos1.lat);
            double deltaLon = Math.abs(geoPos2.lon - geoPos1.lon);
            if (deltaLon > 180) {
                deltaLon += 360;
            }
            Debug.trace("Graticule.create: deltaLat=" + deltaLat + ", deltaLon=" + deltaLon);
            latMajorStep = (float) compose(normalize(gridCellSize * 0.5 * (deltaLon + deltaLat), null));
            lonMajorStep = latMajorStep;
        }
        Debug.trace("Graticule.create: latMajorStep=" + latMajorStep + ", lonMajorStep=" + lonMajorStep);

        float latMinorStep = latMajorStep / 4.0f;
        float lonMinorStep = lonMajorStep / 4.0f;

        int geoBoundaryStep = getGeoBoundaryStep(geoCoding);
        Debug.trace("Graticule.create: geoBoundaryStep=" + geoBoundaryStep);
        final GeoPos[] geoBoundary = ProductUtils.createGeoBoundary(product, null, geoBoundaryStep);
        ProductUtils.normalizeGeoPolygon(geoBoundary);

// nf Debugging, don't delete!
//        GeneralPath generalPath = createPixelBoundaryPath(geoCoding, geoBoundary);
//        if (generalPath != null) {
//            return new Graticule(new GeneralPath[]{generalPath}, null);
//        }

        double xMin = +1.0e10;
        double yMin = +1.0e10;
        double xMax = -1.0e10;
        double yMax = -1.0e10;
        for (int i = 0; i < geoBoundary.length; i++) {
            GeoPos geoPos = geoBoundary[i];
            xMin = Math.min(xMin, geoPos.lon);
            yMin = Math.min(yMin, geoPos.lat);
            xMax = Math.max(xMax, geoPos.lon);
            yMax = Math.max(yMax, geoPos.lat);
        }

        final List parallelList = computeParallelList(product.getGeoCoding(), geoBoundary, latMajorStep, lonMinorStep, yMin, yMax);
        final List meridianList = computeMeridianList(product.getGeoCoding(), geoBoundary, lonMajorStep, latMinorStep, xMin, xMax);
        final GeneralPath[] paths = createPaths(parallelList, meridianList);
        final TextGlyph[] textGlyphs = createTextGlyphs(parallelList, meridianList);

        return new Graticule(paths, textGlyphs);
    }

    /**
     * Creates a graticule for the given product.
     *
     * @param raster               the product
     * @param autoDeterminingSteps if true, <code>gridCellSize</code> is used to compute <code>latMajorStep</code>, <code>lonMajorStep</code> for the given product
     * @param gridCellSize         the grid cell size in pixels, ignored if <code>autoDeterminingSteps</code> if false
     * @param latMajorStep         the grid cell size in meridional direction, ignored if <code>autoDeterminingSteps</code> if true
     * @param lonMajorStep         the grid cell size in parallel direction, ignored if <code>autoDeterminingSteps</code> if true
     * @return the graticule or null, if it could not be created
     */
    public static Graticule create(RasterDataNode raster,
                                   boolean autoDeterminingSteps,
                                   int gridCellSize,
                                   float latMajorStep,
                                   float lonMajorStep) {
        Guardian.assertNotNull("product", raster);
        final GeoCoding geoCoding = raster.getGeoCoding();
        if (geoCoding == null || raster.getSceneRasterWidth() < 16 || raster.getSceneRasterHeight() < 16) {
            return null;
        }

        if (autoDeterminingSteps) {
            final PixelPos pixelPos1 = new PixelPos(0.5f * raster.getSceneRasterWidth(), 0.5f * raster.getSceneRasterHeight());
            final PixelPos pixelPos2 = new PixelPos(pixelPos1.x + 1f, pixelPos1.y + 1f);
            final GeoPos geoPos1 = geoCoding.getGeoPos(pixelPos1, null);
            final GeoPos geoPos2 = geoCoding.getGeoPos(pixelPos2, null);
            double deltaLat = Math.abs(geoPos2.lat - geoPos1.lat);
            double deltaLon = Math.abs(geoPos2.lon - geoPos1.lon);
            if (deltaLon > 180) {
                deltaLon += 360;
            }
            Debug.trace("Graticule.create: deltaLat=" + deltaLat + ", deltaLon=" + deltaLon);
            latMajorStep = (float) compose(normalize(gridCellSize * 0.5 * (deltaLon + deltaLat), null));
            lonMajorStep = latMajorStep;
        }
        Debug.trace("Graticule.create: latMajorStep=" + latMajorStep + ", lonMajorStep=" + lonMajorStep);

        float latMinorStep = latMajorStep / 4.0f;
        float lonMinorStep = lonMajorStep / 4.0f;

        int geoBoundaryStep = getGeoBoundaryStep(geoCoding);
        Debug.trace("Graticule.create: geoBoundaryStep=" + geoBoundaryStep);
        final GeoPos[] geoBoundary = ProductUtils.createGeoBoundary(raster, null, geoBoundaryStep);
        ProductUtils.normalizeGeoPolygon(geoBoundary);

// nf Debugging, don't delete!
//        GeneralPath generalPath = createPixelBoundaryPath(geoCoding, geoBoundary);
//        if (generalPath != null) {
//            return new Graticule(new GeneralPath[]{generalPath}, null);
//        }

        double xMin = +1.0e10;
        double yMin = +1.0e10;
        double xMax = -1.0e10;
        double yMax = -1.0e10;
        for (int i = 0; i < geoBoundary.length; i++) {
            GeoPos geoPos = geoBoundary[i];
            xMin = Math.min(xMin, geoPos.lon);
            yMin = Math.min(yMin, geoPos.lat);
            xMax = Math.max(xMax, geoPos.lon);
            yMax = Math.max(yMax, geoPos.lat);
        }

        final List parallelList = computeParallelList(raster.getGeoCoding(), geoBoundary, latMajorStep, lonMinorStep, yMin, yMax);
        final List meridianList = computeMeridianList(raster.getGeoCoding(), geoBoundary, lonMajorStep, latMinorStep, xMin, xMax);
        final GeneralPath[] paths = createPaths(parallelList, meridianList);
        final TextGlyph[] textGlyphs = createTextGlyphs(parallelList, meridianList);

        return new Graticule(paths, textGlyphs);
    }

    private static int getGeoBoundaryStep(final GeoCoding geoCoding) {
        int step = 16;
        if (geoCoding instanceof TiePointGeoCoding) {
            final TiePointGeoCoding tiePointGeoCoding = (TiePointGeoCoding) geoCoding;
            step = Math.round(Math.min(tiePointGeoCoding.getLonGrid().getSubSamplingX(), tiePointGeoCoding.getLonGrid().getSubSamplingY()));
        }
        return step;
    }

    private static ArrayList computeParallelList(final GeoCoding geoCoding, final GeoPos[] geoBoundary,
                                                 final double latMajorStep,
                                                 final double lonMinorStep,
                                                 final double yMin,
                                                 final double yMax) {
//        final GeoCoding geoCoding = product.getGeoCoding();
        ArrayList parallelList = new ArrayList();
        ArrayList intersectionList = new ArrayList();
        GeoPos geoPos, int1, int2;
        PixelPos pixelPos;
        float lat, lon;
        double my = latMajorStep * Math.floor(yMin / latMajorStep);
        for (; my <= yMax; my += latMajorStep) {
            intersectionList.clear();
            computeParallelIntersections(geoBoundary, my, intersectionList);
            if (intersectionList.size() > 0 && intersectionList.size() % 2 == 0) {
                final GeoPos[] intersections = (GeoPos[]) intersectionList.toArray(new GeoPos[intersectionList.size()]);
                Arrays.sort(intersections, new GeoPosLonComparator());
                ArrayList parallel = new ArrayList();
                // loop forward order
                for (int i = 0; i < intersections.length; i += 2) {
                    int1 = intersections[i];
                    int2 = intersections[i + 1];
                    lat = int1.lat;
                    lon = int1.lon;
                    for (int k = 0; k <= 1;) {
                        geoPos = new GeoPos(lat, limitLon(lon));
                        pixelPos = geoCoding.getPixelPos(geoPos, null);
                        parallel.add(new Coord(geoPos, pixelPos));
                        lon += lonMinorStep;
                        if (lon >= int2.lon) {
                            lon = int2.lon;
                            k++;
                        }
                    }
                }
                parallelList.add(parallel);
            }
        }
        return parallelList;
    }

    private static ArrayList computeMeridianList(final GeoCoding geoCoding, final GeoPos[] geoBoundary,
                                                 final double lonMajorStep,
                                                 final double latMinorStep,
                                                 final double xMin,
                                                 final double xMax) {
//        final GeoCoding geoCoding = product.getGeoCoding();
        ArrayList meridianList = new ArrayList();
        ArrayList intersectionList = new ArrayList();
        GeoPos geoPos, int1, int2;
        PixelPos pixelPos;
        float lat, lon;
        double mx = lonMajorStep * Math.floor(xMin / lonMajorStep);
        for (; mx <= xMax; mx += lonMajorStep) {
            intersectionList.clear();
            computeMeridianIntersections(geoBoundary, mx, intersectionList);
            if (intersectionList.size() > 0 && intersectionList.size() % 2 == 0) {
                final GeoPos[] intersections = (GeoPos[]) intersectionList.toArray(new GeoPos[intersectionList.size()]);
                Arrays.sort(intersections, new GeoPosLatComparator());
                ArrayList meridian = new ArrayList();
                // loop reverse order
                for (int i = intersections.length - 2; i >= 0; i -= 2) {
                    int1 = intersections[i + 1];
                    int2 = intersections[i];
                    lat = int1.lat;
                    lon = int1.lon;
                    for (int k = 0; k <= 1;) {
                        geoPos = new GeoPos(lat, limitLon(lon));
                        pixelPos = geoCoding.getPixelPos(geoPos, null);
                        meridian.add(new Coord(geoPos, pixelPos));
                        lat -= latMinorStep;
                        if (lat <= int2.lat) {
                            lat = int2.lat;
                            k++;
                        }
                    }
                }
                meridianList.add(meridian);
            }
        }
        return meridianList;
    }

    private static void computeParallelIntersections(final GeoPos[] geoBoundary,
                                                     final double my,
                                                     final List intersectionList) {
        double p0x = 0, p0y = 0;
        double p1x, p1y;
        double pa;
        double mx;
        for (int i = 0; i < geoBoundary.length; i++) {
            GeoPos geoPos = geoBoundary[i];
            p1x = geoPos.lon;
            p1y = geoPos.lat;
            if (i > 0) {
                if (((my >= p0y && my <= p1y) || (my >= p1y && my <= p0y)) &&
                        (p1y - p0y != 0.0)) {
                    pa = (my - p0y) / (p1y - p0y);
                    if (pa >= 0.0 && pa < 1.0) {
                        mx = p0x + pa * (p1x - p0x);
                        intersectionList.add(new GeoPos((float) my, (float) mx));
                    }
                }
            }
            p0x = p1x;
            p0y = p1y;
        }
    }

    private static void computeMeridianIntersections(final GeoPos[] geoBoundary,
                                                     final double mx,
                                                     final List intersectionList) {
        double p0x = 0, p0y = 0;
        double p1x, p1y;
        double pa;
        double my;
        for (int i = 0; i < geoBoundary.length; i++) {
            GeoPos geoPos = geoBoundary[i];
            p1x = geoPos.lon;
            p1y = geoPos.lat;
            if (i > 0) {
                if (((mx >= p0x && mx <= p1x) || (mx >= p1x && mx <= p0x)) &&
                        (p1x - p0x != 0.0)) {
                    pa = (mx - p0x) / (p1x - p0x);
                    if (pa >= 0.0 && pa < 1.0) {
                        my = p0y + pa * (p1y - p0y);
                        intersectionList.add(new GeoPos((float) my, (float) mx));
                    }
                }
            }
            p0x = p1x;
            p0y = p1y;
        }
    }

    private static GeneralPath[] createPaths(List parallelList, List meridianList) {
        ArrayList generalPathList = new ArrayList();
        addToPath(parallelList, generalPathList);
        addToPath(meridianList, generalPathList);
        return (GeneralPath[]) generalPathList.toArray(new GeneralPath[generalPathList.size()]);
    }

    private static void addToPath(List lineList, List generalPathList) {
        for (int i = 0; i < lineList.size(); i++) {
            List coordList = (List) lineList.get(i);
            if (coordList.size() >= 2) {
                final GeneralPath generalPath = new GeneralPath();
                for (int j = 0; j < coordList.size(); j++) {
                    PixelPos pixelPos = ((Coord) coordList.get(j)).pixelPos;
                    if (j == 0) {
                        generalPath.moveTo(pixelPos.x, pixelPos.y);
                    } else {
                        generalPath.lineTo(pixelPos.x, pixelPos.y);
                    }
                }
                generalPathList.add(generalPath);
            }
        }
    }

    private static TextGlyph[] createTextGlyphs(List parallelList, List meridianList) {
        ArrayList textGlyphList = new ArrayList();
        createParallelTextGlyphs(parallelList, textGlyphList);
        createMeridianTextGlyphs(meridianList, textGlyphList);
        return (TextGlyph[]) textGlyphList.toArray(new TextGlyph[textGlyphList.size()]);
    }


    private static void createParallelTextGlyphs(List parallelList,
                                                 List textGlyphList) {
        Coord coord1;
        Coord coord2;
        for (int i = 0; i < parallelList.size(); i++) {
            final List parallel = (List) parallelList.get(i);
            if (parallel.size() >= 3) {
                coord1 = (Coord) parallel.get(1);
                coord2 = (Coord) parallel.get(2);
                textGlyphList.add(createLatTextGlyph(coord1, coord2));
            } else if (parallel.size() >= 2) {
                coord1 = (Coord) parallel.get(0);
                coord2 = (Coord) parallel.get(1);
                textGlyphList.add(createLatTextGlyph(coord1, coord2));
            }
        }
    }

    private static void createMeridianTextGlyphs(List meridianList,
                                                 List textGlyphList) {
        Coord coord1;
        Coord coord2;
        for (int i = 0; i < meridianList.size(); i++) {
            final List meridian = (List) meridianList.get(i);
            if (meridian.size() >= 3) {
                coord1 = (Coord) meridian.get(1);
                coord2 = (Coord) meridian.get(2);
                textGlyphList.add(createLonTextGlyph(coord1, coord2));
            } else if (meridian.size() >= 2) {
                coord1 = (Coord) meridian.get(0);
                coord2 = (Coord) meridian.get(1);
                textGlyphList.add(createLonTextGlyph(coord1, coord2));
            }
        }
    }

    private static TextGlyph createLatTextGlyph(Coord coord1, Coord coord2) {
        return createTextGlyph(coord1.geoPos.getLatString(), coord1, coord2);
    }

    private static TextGlyph createLonTextGlyph(Coord coord1, Coord coord2) {
        return createTextGlyph(coord1.geoPos.getLonString(), coord1, coord2);
    }

    private static TextGlyph createTextGlyph(String text, Coord coord1, Coord coord2) {
        final float angle = (float) Math.atan2(coord2.pixelPos.y - coord1.pixelPos.y,
                                               coord2.pixelPos.x - coord1.pixelPos.x);
        return new TextGlyph(text, coord1.pixelPos.x, coord1.pixelPos.y, angle);
    }

    private static float limitLon(float lon) {
        while (lon < -180f) {
            lon += 360f;
        }
        while (lon > 180f) {
            lon -= 360f;
        }
        return lon;
    }

    private static double[] normalize(double x, double[] result) {
        final double exponent = (x == 0.0) ? 0.0 : Math.ceil(Math.log(Math.abs(x)) / Math.log(10.0));
        final double mantissa = (x == 0.0) ? 0.0 : x / Math.pow(10.0, exponent);
        if (result == null) {
            result = new double[2];
        }
        result[0] = mantissa;
        result[1] = exponent;
        return result;
    }

    private static double compose(final double[] components) {
        final double mantissa = components[0];
        final double exponent = components[1];
        final double mantissaRounded;
        if (mantissa < 0.15) {
            mantissaRounded = 0.1;
        } else if (mantissa < 0.225) {
            mantissaRounded = 0.2;
        } else if (mantissa < 0.375) {
            mantissaRounded = 0.25;
        } else if (mantissa < 0.75) {
            mantissaRounded = 0.5;
        } else {
            mantissaRounded = 1.0;
        }
        return mantissaRounded * Math.pow(10.0, exponent);
    }

    /**
     * Not used, but usefull for debugging: DON'T delete this method!
     *
     * @param geoCoding
     * @param geoBoundary
     * @return the geo-boundary
     */
    private static GeneralPath createPixelBoundaryPath(final GeoCoding geoCoding, final GeoPos[] geoBoundary) {
        final GeneralPath generalPath = new GeneralPath();
        for (int i = 0; i < geoBoundary.length; i++) {
            final GeoPos geoPos = geoBoundary[i];
            geoPos.lon = limitLon(geoPos.lon);
            final PixelPos pixelPos = geoCoding.getPixelPos(geoPos, null);
            if (i == 0) {
                generalPath.moveTo(pixelPos.x, pixelPos.y);
            } else {
                generalPath.lineTo(pixelPos.x, pixelPos.y);
            }
        }
        return generalPath;
    }

    public static class TextGlyph {

        private final String text;
        private final float x;
        private final float y;
        private final float angle;

        public TextGlyph(String text, float x, float y, float angle) {
            this.text = text;
            this.x = x;
            this.y = y;
            this.angle = angle;
        }

        public String getText() {
            return text;
        }

        public float getX() {
            return x;
        }

        public float getY() {
            return y;
        }

        public float getAngle() {
            return angle;
        }
    }

    private static class Coord {
        GeoPos geoPos;
        PixelPos pixelPos;

        public Coord(GeoPos geoPos, PixelPos pixelPos) {
            this.geoPos = geoPos;
            this.pixelPos = pixelPos;
        }
    }

    private static class GeoPosLatComparator implements Comparator {
        public int compare(Object o1, Object o2) {
            final GeoPos geoPos1 = (GeoPos) o1;
            final GeoPos geoPos2 = (GeoPos) o2;
            final float delta = geoPos1.lat - geoPos2.lat;
            if (delta < 0f) {
                return -1;
            } else if (delta > 0f) {
                return 1;
            } else {
                return 0;
            }
        }
    }

    private static class GeoPosLonComparator implements Comparator {
        public int compare(Object o1, Object o2) {
            final GeoPos geoPos1 = (GeoPos) o1;
            final GeoPos geoPos2 = (GeoPos) o2;
            final float delta = geoPos1.lon - geoPos2.lon;
            if (delta < 0f) {
                return -1;
            } else if (delta > 0f) {
                return 1;
            } else {
                return 0;
            }
        }
    }
}
