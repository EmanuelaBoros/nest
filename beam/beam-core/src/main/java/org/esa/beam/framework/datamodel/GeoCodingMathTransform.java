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
 * @version $Revision: 1.1 $ $Date: 2009-04-28 14:39:32 $
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

        if (!geoCoding.equals(that.geoCoding)) {
            return false;
        }
        return t == that.t;

    }

    @Override
    public int hashCode() {
        int result = super.hashCode();
        result = 31 * result + geoCoding.hashCode();
        result = 31 * result + t.hashCode();
        return result;
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
        }
    }

    private static class TG2P implements T {

        @Override
        public void transform(GeoCoding geoCoding,
                              double[] srcPts, int srcOff,
                              double[] dstPts, int dstOff,
                              int numPts) throws TransformException {
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
        }
    }
}