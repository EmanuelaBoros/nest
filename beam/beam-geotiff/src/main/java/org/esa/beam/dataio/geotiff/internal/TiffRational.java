package org.esa.beam.dataio.geotiff.internal;

import org.esa.beam.framework.datamodel.ProductData;

/**
 * A TIFFValue implementation for the GeoTIFF format.
 *
 * @author Marco Peters
 * @author Sabine Embacher
 * @author Norman Fomferra
 * @version $Revision: 1.1 $ $Date: 2009-04-28 14:37:14 $
 */
class TiffRational extends TiffValue {

    private static final int NUMERATOR_INDEX = 0;
    private static final int DENOMINATOR_INDEX = 1;

    public TiffRational(final long numerator, final long denominator) {
        TiffValueRangeChecker.checkValueTiffRational(numerator, "numerator");
        TiffValueRangeChecker.checkValueTiffRational(denominator, "denominator");
        setData(ProductData.createInstance(ProductData.TYPE_UINT32, 2));
        getData().setElemUIntAt(NUMERATOR_INDEX, numerator);
        getData().setElemUIntAt(DENOMINATOR_INDEX, denominator);
    }

    public long getNumerator() {
        return getData().getElemUIntAt(NUMERATOR_INDEX);
    }

    public long getDenominator() {
        return getData().getElemUIntAt(DENOMINATOR_INDEX);
    }

    public double getValue() {
        return getNumerator() / getDenominator();
    }
}
