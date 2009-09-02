package org.esa.beam.framework.datamodel;

import org.geotools.parameter.DefaultParameterDescriptorGroup;
import org.geotools.referencing.operation.transform.AbstractMathTransform;
import org.opengis.parameter.GeneralParameterDescriptor;
import org.opengis.parameter.ParameterDescriptorGroup;
import org.opengis.referencing.operation.MathTransform;
import org.opengis.referencing.operation.NoninvertibleTransformException;
import org.opengis.referencing.operation.TransformException;

/**
 * A math transform which converts from grid (pixel) coordinates to geograhical coordinates.
 *
 * @author Marco Peters
 * @author Norman Fomferra
 * @version $Revision: 1.3 $ $Date: 2009-09-01 20:27:12 $
 * @since BEAM 4.6
 */
public class GeoCodingMathTransform extends AbstractMathTransform {

    private static final TG2P G2P = new TG2P();
    private static final TP2G P2G = new TP2G();
    private static final int DIMS = 2;

    private final GeoCoding geoCoding;
    private final T t;


    public GeoCodingMathTransform(GeoCoding geoCoding) {
        this(geoCoding, G2P);
    }

    private GeoCodingMathTransform(GeoCoding geoCoding, T t) {
        this.geoCoding = geoCoding;
        this.t = t;
    }

    @Override
    public ParameterDescriptorGroup getParameterDescriptors() {
        return new DefaultParameterDescriptorGroup(getClass().getSimpleName(), new GeneralParameterDescriptor[0]);

    }

    @Override
    public int getSourceDimensions() {
        return DIMS;
    }

    @Override
    public int getTargetDimensions() {
        return DIMS;
    }

    @Override
    public MathTransform inverse() throws NoninvertibleTransformException {
        return new GeoCodingMathTransform(geoCoding, t == G2P ? P2G : G2P);
    }

    @Override
    public void transform(double[] srcPts, int srcOff,
                          double[] dstPts, int dstOff,
                          int numPts) throws TransformException {
        t.transform(geoCoding, srcPts, srcOff, dstPts, dstOff, numPts);
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (o == null || getClass() != o.getClass()) {
            return false;
        }
        if (!super.equals(o)) {
            return false;
        }
        
        GeoCodingMathTransform that = (GeoCodingMathTransform) o;
        if (t != that.t) {
            return false;
        }

        if (areGeoCodingsEqual(geoCoding, that.geoCoding)) {
            return true;
        }
        
        return false;
    }

    @Override
    public int hashCode() {
        int result = super.hashCode();
        result = 31 * result + geoCoding.hashCode();
        result = 31 * result + t.hashCode();
        return result;
    }
    
    private static boolean areGeoCodingsEqual(GeoCoding oneGeoCoding, GeoCoding otherGeoCoding) {
        try {
            final PixelPos pixelPos0 = new PixelPos(0, 0);
            final GeoPos geoPos = oneGeoCoding.getGeoPos(pixelPos0, null);
            final GeoPos thatGeoPos = otherGeoCoding.getGeoPos(pixelPos0, null);
            if (!equalsGeoPos(geoPos, thatGeoPos)) {
                return false;
            }
            final PixelPos pixelPos1 = new PixelPos(1, 1);
            final GeoPos geoPos1 = oneGeoCoding.getGeoPos(pixelPos1, null);
            final GeoPos thatGeoPos1 = otherGeoCoding.getGeoPos(pixelPos1, null);
            return equalsGeoPos(geoPos1, thatGeoPos1);
        } catch (Exception e) {
            return false;
        }
    }

    private static boolean equalsGeoPos(GeoPos geoPos, GeoPos thatGeoPos) {
        return Math.abs(geoPos.getLat() - thatGeoPos.getLat()) < 1.0e-6 &&
               Math.abs(geoPos.getLon() - thatGeoPos.getLon()) < 1.0e-6;
    }

    private interface T {

        void transform(GeoCoding geoCoding, double[] srcPts, int srcOff,
                       double[] dstPts, int dstOff,
                       int numPts) throws TransformException;
    }

    private static class TP2G implements T {

        @Override
        public void transform(GeoCoding geoCoding,
                              double[] srcPts, int srcOff,
                              double[] dstPts, int dstOff,
                              int numPts) throws TransformException {
            try {
                GeoPos geoPos = new GeoPos();
                PixelPos pixelPos = new PixelPos();
                for (int i = 0; i < numPts; i++) {
                    final int firstIndex = (DIMS * i);
                    final int secondIndex = firstIndex + 1;
                    pixelPos.x = (float) srcPts[srcOff + firstIndex];
                    pixelPos.y = (float) srcPts[srcOff + secondIndex];

                    geoCoding.getGeoPos(pixelPos, geoPos);

                    dstPts[dstOff + firstIndex] = geoPos.lon;
                    dstPts[dstOff + secondIndex] = geoPos.lat;
                }
            } catch (Exception e) {
                TransformException transformException = new TransformException();
                transformException.initCause(e);
                throw transformException;
            }
        }
    }

    private static class TG2P implements T {

        @Override
        public void transform(GeoCoding geoCoding,
                              double[] srcPts, int srcOff,
                              double[] dstPts, int dstOff,
                              int numPts) throws TransformException {
            try {
                GeoPos geoPos = new GeoPos();
                PixelPos pixelPos = new PixelPos();
                for (int i = 0; i < numPts; i++) {
                    final int firstIndex = (DIMS * i);
                    final int secondIndex = firstIndex + 1;
                    geoPos.lon = (float) srcPts[srcOff + firstIndex];
                    geoPos.lat = (float) srcPts[srcOff + secondIndex];

                    geoCoding.getPixelPos(geoPos, pixelPos);

                    dstPts[dstOff + firstIndex] = pixelPos.x;
                    dstPts[dstOff + secondIndex] = pixelPos.y;
                }
            } catch (Exception e) {
                TransformException transformException = new TransformException();
                transformException.initCause(e);
                throw transformException;
            }
        }
    }
}
