/*
 * $Id: GeneralFilterBand.java,v 1.3 2009-05-28 14:17:58 lveci Exp $
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
import com.bc.ceres.glevel.MultiLevelModel;
import com.bc.ceres.glevel.support.AbstractMultiLevelSource;
import com.bc.ceres.glevel.support.DefaultMultiLevelImage;
import org.esa.beam.jai.ImageManager;

import javax.media.jai.BorderExtender;
import javax.media.jai.BorderExtenderCopy;
import javax.media.jai.JAI;
import javax.media.jai.KernelJAI;
import javax.media.jai.operator.ConvolveDescriptor;
import javax.media.jai.operator.MaxFilterDescriptor;
import javax.media.jai.operator.MaxFilterShape;
import javax.media.jai.operator.MedianFilterDescriptor;
import javax.media.jai.operator.MedianFilterShape;
import javax.media.jai.operator.MinFilterDescriptor;
import javax.media.jai.operator.MinFilterShape;
import java.awt.RenderingHints;
import java.awt.image.RenderedImage;
import java.io.IOException;
import java.util.Arrays;

/**
 * A band that obtains its input data from an underlying source raster and filters
 * its data using an arbitrary {@link Operator algorithm}.
 * <p/>
 * <p><i>Note that this class is not yet public API and may change in future releases.</i></p>
 *
 * @author Norman Fomferra
 * @version $Revision: 1.3 $ $Date: 2009-05-28 14:17:58 $
 */
public class GeneralFilterBand extends FilterBand {

    public static final Operator MIN = new Min();               // JAI: MinFilterDescriptor
    public static final Operator MAX = new Max();               // JAI: MaxFilterDescriptor
    public static final Operator MEDIAN = new Median();         // JAI: MedianFilterDescriptor
    public static final Operator MEAN = new Mean();              // JAI: ConvolveDescriptor
    public static final Operator STDDEV = new StandardDeviation();     // TODO - Write JAI Operator
    public static final Operator RMS = new RootMeanSquare();           // TODO - Write JAI Operator

    private final int subWindowSize;
    private final Operator operator;

    /**
     * Creates a GeneralFilterBand.
     *
     * @param name            the name of the band.
     * @param source          the source which shall be filtered.
     * @param subWindowWidth  the window width used by the filter
     * @param subWindowHeight the window height used by the filter
     * @param operator        the operator which performs the filter operation
     * @deprecated since BEAM 4.5; non square windows are not supported, use {@link GeneralFilterBand#GeneralFilterBand(String, RasterDataNode, int, Operator)}.
     *             It is currently implemented by delegating to the other constructor using {@code subWindowWidth}  as {@code subWindowSize}.
     */
    @Deprecated
    public GeneralFilterBand(String name, RasterDataNode source, int subWindowWidth, int subWindowHeight,
                             Operator operator) {
        this(name, source, subWindowWidth, operator);
    }

    /**
     * Creates a GeneralFilterBand.
     *
     * @param name          the name of the band.
     * @param source        the source which shall be filtered.
     * @param subWindowSize the window size (width/height) used by the filter
     * @param operator      the operator which performs the filter operation
     */
    public GeneralFilterBand(String name, RasterDataNode source, int subWindowSize, Operator operator) {
        super(name,
              source.getGeophysicalDataType(),
              source.getSceneRasterWidth(),
              source.getSceneRasterHeight(),
              source);
        setOwner(source.getProduct());
        this.subWindowSize = subWindowSize;
        this.operator = operator;
    }

    @Override
    protected RenderedImage createSourceImage() {
        final MultiLevelModel multiLevelModel = ImageManager.getInstance().getMultiLevelModel(this);
        return new DefaultMultiLevelImage(new GeneralFilterMultiLevelSource(multiLevelModel, getSource(),
                                                                            BorderExtender.createInstance(
                                                                                    BorderExtenderCopy.BORDER_COPY)));
    }

    /**
     * Creates an instance {@link Operator} by the given class name
     *
     * @param operatorClassName the class name
     * @return instance of {@link Operator}
     */
    public static Operator createOperator(String operatorClassName) {
        Operator operator = null;
        try {
            final Class operatorClass = Class.forName(operatorClassName);
            operator = (Operator) operatorClass.newInstance();
        } catch (ClassNotFoundException e) {
            e.printStackTrace();
        } catch (IllegalAccessException e) {
            e.printStackTrace();
        } catch (InstantiationException e) {
            e.printStackTrace();
        }
        return operator;
    }

    public int getSubWindowSize() {
        return subWindowSize;
    }

    /**
     * @return return the width of the sub window.
     * @deprecated since BEAM 4.5, use {@link #getSubWindowSize()}
     */
    @Deprecated
    public int getSubWindowWidth() {
        return subWindowSize;
    }

    /**
     * @return return the height of the sub window.
     * @deprecated since BEAM 4.5, use {@link #getSubWindowSize()}
     */
    @Deprecated
    public int getSubWindowHeight() {
        return subWindowSize;
    }

    public Operator getOperator() {
        return operator;
    }

    /**
     * An operator which performs an operation on an array of pixel values extracted from a raster sub-window.
     */
    public static interface Operator {

        /**
         * Performs the operation.
         *
         * @param values an array of pixel values extracted from a raster sub-window
         * @param n      the number of values in the array
         * @return the result of the operation
         */
        double evaluate(double[] values, int n);
    }

    public static class Min implements Operator {

        public final double evaluate(final double[] values, final int n) {
            if (n == 0) {
                return Double.NaN;
            }
            double x = Double.MAX_VALUE;
            for (int i = 0; i < n; i++) {
                if (values[i] < x) {
                    x = values[i];
                }
            }
            return x;
        }
    }

    public static class Max implements Operator {

        public final double evaluate(final double[] values, final int n) {
            if (n == 0) {
                return Double.NaN;
            }
            double x = -Double.MAX_VALUE;
            for (int i = 0; i < n; i++) {
                if (values[i] > x) {
                    x = values[i];
                }
            }
            return x;
        }
    }

    public static class Mean implements Operator {

        public final double evaluate(final double[] values, final int n) {
            if (n == 0) {
                return Double.NaN;
            }
            double sum = 0;
            for (int i = 0; i < n; i++) {
                sum += values[i];
            }
            return sum / n;
        }
    }

    public static class Median implements Operator {

        public final double evaluate(final double[] values, final int n) {
            if (n == 0) {
                return Double.NaN;
            } else if (n == 1) {
                return values[0];
            }
            Arrays.sort(values);
            final int i = n / 2;
            if (n % 2 != 0) {
                return values[i];
            }
            return 0.5 * (values[i] + values[i + 1]);
        }
    }

    public static class StandardDeviation implements Operator {

        public final double evaluate(final double[] values, final int n) {
            if (n == 0) {
                return Double.NaN;
            } else if (n == 1) {
                return 0;
            }
            double sum = 0;
            for (int i = 0; i < n; i++) {
                sum += values[i];
            }
            final double mean = sum / n;
            sum = 0;
            for (int i = 0; i < n; i++) {
                sum += (values[i] - mean) * (values[i] - mean);
            }
            return Math.sqrt(sum / (n - 1));
        }
    }

    public static class RootMeanSquare implements Operator {

        public final double evaluate(final double[] values, final int n) {
            if (n == 0) {
                return Double.NaN;
            }
            double sum = 0;
            for (int i = 0; i < n; i++) {
                sum += values[i] * values[i];
            }
            return Math.sqrt(sum / n);
        }
    }

    private class GeneralFilterMultiLevelSource extends AbstractMultiLevelSource {

        private final RasterDataNode sourceRaster;
        private BorderExtender noDataExtender;

        private GeneralFilterMultiLevelSource(final MultiLevelModel multiLevelModel, RasterDataNode sourceRaster,
                                              final BorderExtender borderExtender) {
            super(multiLevelModel);
            this.sourceRaster = sourceRaster;
            noDataExtender = borderExtender;
        }

        @Override
        public RenderedImage createImage(int level) {
            final ImageManager imageManager = ImageManager.getInstance();
            final RenderedImage geophysicalImage = imageManager.getGeophysicalImage(sourceRaster, level);
            if (getOperator() == MIN) {
                final MinFilterShape maskSquare = MinFilterDescriptor.MIN_MASK_SQUARE;
                return MinFilterDescriptor.create(geophysicalImage, maskSquare, subWindowSize, null);
            }
            if (getOperator() == MAX) {
                final MaxFilterShape maskSquare = MaxFilterDescriptor.MAX_MASK_SQUARE;
                return MaxFilterDescriptor.create(geophysicalImage, maskSquare, subWindowSize, null);
            }
            if (getOperator() == MEDIAN) {
                final MedianFilterShape maskSquare = MedianFilterDescriptor.MEDIAN_MASK_SQUARE;
                return MedianFilterDescriptor.create(geophysicalImage, maskSquare, subWindowSize, null);
            }
            if (getOperator() == MEAN) {
                final int kernelSize = subWindowSize * subWindowSize;
                float[] meanFilter = new float[kernelSize];
                Arrays.fill(meanFilter, 1.0f / kernelSize);
                int keyOrigin = (int) Math.ceil(subWindowSize / 2.0f);
                KernelJAI kernel = new KernelJAI(subWindowSize, subWindowSize, keyOrigin, keyOrigin, meanFilter);
                RenderingHints rh = new RenderingHints(JAI.KEY_BORDER_EXTENDER, noDataExtender);
                return ConvolveDescriptor.create(geophysicalImage, kernel, rh);
            }
            if (getOperator() == STDDEV) {
                // TODO - implement a JAI operator
            }
            if (getOperator() == RMS) {
                // TODO - implement a JAI operator
            }
            throw new IllegalStateException(
                    String.format("Operator class %s not supported.", getOperator().getClass()));

        }

    }
}