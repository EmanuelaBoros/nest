package org.esa.beam.framework.datamodel;

import javax.media.jai.PixelAccessor;
import javax.media.jai.UnpackedImageData;
import java.awt.image.Raster;
import java.awt.image.DataBuffer;
import java.awt.Rectangle;

/**
*
* @author Norman Fomferra
* @author Marco Peters
* @version $Revision: 1.2 $ $Date: 2010-02-12 22:18:10 $
* @since BEAM 4.5.1
*/
class HistogramStxOp implements StxOp {

    private final double lowValue;
    private final double highValue;
    private final double binWidth;

    private final int[] bins;

    HistogramStxOp(int numBins, double lowValue, double highValue) {
        this.lowValue = lowValue;
        this.highValue = highValue;
        this.binWidth = (highValue - lowValue) / numBins;
        this.bins = new int[numBins];
    }

    public String getName() {
        return "Histogram";
    }

    public int[] getBins() {
        return bins;
    }

    public void accumulateDataUByte(PixelAccessor dataAccessor, Raster dataTile, PixelAccessor maskAccessor,
                                    Raster maskTile, Rectangle r, String unit) {
        final int[] bins = this.bins;
        final double lowValue = this.lowValue;
        final double highValue = this.highValue;
        final double binWidth = this.binWidth;

        final UnpackedImageData duid = dataAccessor.getPixels(dataTile, r, DataBuffer.TYPE_BYTE, false);
        final byte[] data = duid.getByteData(0);
        final int dataPixelStride = duid.pixelStride;
        final int dataLineStride = duid.lineStride;
        final int dataBandOffset = duid.bandOffsets[0];

        byte[] mask = null;
        int maskPixelStride = 0;
        int maskLineStride = 0;
        int maskBandOffset = 0;
        if (maskAccessor != null) {
            UnpackedImageData muid = maskAccessor.getPixels(maskTile, r, DataBuffer.TYPE_BYTE, false);
            mask = muid.getByteData(0);
            maskPixelStride = muid.pixelStride;
            maskLineStride = muid.lineStride;
            maskBandOffset = muid.bandOffsets[0];
        }

        final int width = r.width;
        final int height = r.height;

        int dataLineOffset = dataBandOffset;
        int maskLineOffset = maskBandOffset;
        for (int y = 0; y < height; y++) {
            int dataPixelOffset = dataLineOffset;
            int maskPixelOffset = maskLineOffset;
            for (int x = 0; x < width; x++) {
                if (mask == null || mask[maskPixelOffset] != 0) {
                    double d = data[dataPixelOffset] & 0xff;
                    if (d >= lowValue && d < highValue) {
                        int i = (int) ((d - lowValue) / binWidth);
                        bins[i]++;
                    }
                }
                dataPixelOffset += dataPixelStride;
                maskPixelOffset += maskPixelStride;
            }
            dataLineOffset += dataLineStride;
            maskLineOffset += maskLineStride;
        }
    }

    public void accumulateDataUShort(PixelAccessor dataAccessor, Raster dataTile, PixelAccessor maskAccessor,
                                     Raster maskTile, Rectangle r, String unit) {
        final int[] bins = this.bins;
        final double lowValue = this.lowValue;
        final double highValue = this.highValue;
        final double binWidth = this.binWidth;

        final UnpackedImageData duid = dataAccessor.getPixels(dataTile, r, DataBuffer.TYPE_USHORT, false);
        final short[] data = duid.getShortData(0);
        final int dataPixelStride = duid.pixelStride;
        final int dataLineStride = duid.lineStride;
        final int dataBandOffset = duid.bandOffsets[0];

        byte[] mask = null;
        int maskPixelStride = 0;
        int maskLineStride = 0;
        int maskBandOffset = 0;
        if (maskAccessor != null) {
            UnpackedImageData muid = maskAccessor.getPixels(maskTile, r, DataBuffer.TYPE_BYTE, false);
            mask = muid.getByteData(0);
            maskPixelStride = muid.pixelStride;
            maskLineStride = muid.lineStride;
            maskBandOffset = muid.bandOffsets[0];
        }

        final int width = r.width;
        final int height = r.height;

        int dataLineOffset = dataBandOffset;
        int maskLineOffset = maskBandOffset;
        for (int y = 0; y < height; y++) {
            int dataPixelOffset = dataLineOffset;
            int maskPixelOffset = maskLineOffset;
            for (int x = 0; x < width; x++) {
                if (mask == null || mask[maskPixelOffset] != 0) {
                    double d = data[dataPixelOffset] & 0xffff;
                    if (d >= lowValue && d < highValue) {
                        int i = (int) ((d - lowValue) / binWidth);
                        bins[i]++;
                    }
                }
                dataPixelOffset += dataPixelStride;
                maskPixelOffset += maskPixelStride;
            }
            dataLineOffset += dataLineStride;
            maskLineOffset += maskLineStride;
        }
    }

    public void accumulateDataShort(PixelAccessor dataAccessor, Raster dataTile, PixelAccessor maskAccessor,
                                    Raster maskTile, Rectangle r, String unit) {
        final int[] bins = this.bins;
        final double lowValue = this.lowValue;
        final double highValue = this.highValue;
        final double binWidth = this.binWidth;

        final UnpackedImageData duid = dataAccessor.getPixels(dataTile, r, DataBuffer.TYPE_SHORT, false);
        final short[] data = duid.getShortData(0);
        final int dataPixelStride = duid.pixelStride;
        final int dataLineStride = duid.lineStride;
        final int dataBandOffset = duid.bandOffsets[0];

        byte[] mask = null;
        int maskPixelStride = 0;
        int maskLineStride = 0;
        int maskBandOffset = 0;
        if (maskAccessor != null) {
            UnpackedImageData muid = maskAccessor.getPixels(maskTile, r, DataBuffer.TYPE_BYTE, false);
            mask = muid.getByteData(0);
            maskPixelStride = muid.pixelStride;
            maskLineStride = muid.lineStride;
            maskBandOffset = muid.bandOffsets[0];
        }

        final int width = r.width;
        final int height = r.height;

        int dataLineOffset = dataBandOffset;
        int maskLineOffset = maskBandOffset;
        for (int y = 0; y < height; y++) {
            int dataPixelOffset = dataLineOffset;
            int maskPixelOffset = maskLineOffset;
            for (int x = 0; x < width; x++) {
                if (mask == null || mask[maskPixelOffset] != 0) {
                    double d = data[dataPixelOffset];
                    if (d >= lowValue && d < highValue) {
                        int i = (int) ((d - lowValue) / binWidth);
                        bins[i]++;
                    }
                }
                dataPixelOffset += dataPixelStride;
                maskPixelOffset += maskPixelStride;
            }
            dataLineOffset += dataLineStride;
            maskLineOffset += maskLineStride;
        }
    }

    public void accumulateDataInt(PixelAccessor dataAccessor, Raster dataTile, PixelAccessor maskAccessor,
                                  Raster maskTile, Rectangle r, String unit) {
        final int[] bins = this.bins;
        final double lowValue = this.lowValue;
        final double highValue = this.highValue;
        final double binWidth = this.binWidth;

        final UnpackedImageData duid = dataAccessor.getPixels(dataTile, r, DataBuffer.TYPE_INT, false);
        final int[] data = duid.getIntData(0);
        final int dataPixelStride = duid.pixelStride;
        final int dataLineStride = duid.lineStride;
        final int dataBandOffset = duid.bandOffsets[0];

        byte[] mask = null;
        int maskPixelStride = 0;
        int maskLineStride = 0;
        int maskBandOffset = 0;
        if (maskAccessor != null) {
            UnpackedImageData muid = maskAccessor.getPixels(maskTile, r, DataBuffer.TYPE_BYTE, false);
            mask = muid.getByteData(0);
            maskPixelStride = muid.pixelStride;
            maskLineStride = muid.lineStride;
            maskBandOffset = muid.bandOffsets[0];
        }

        final int width = r.width;
        final int height = r.height;

        int dataLineOffset = dataBandOffset;
        int maskLineOffset = maskBandOffset;
        for (int y = 0; y < height; y++) {
            int dataPixelOffset = dataLineOffset;
            int maskPixelOffset = maskLineOffset;
            for (int x = 0; x < width; x++) {
                if (mask == null || mask[maskPixelOffset] != 0) {
                    double d = data[dataPixelOffset];
                    if (d >= lowValue && d < highValue) {
                        int i = (int) ((d - lowValue) / binWidth);
                        bins[i]++;
                    }
                }
                dataPixelOffset += dataPixelStride;
                maskPixelOffset += maskPixelStride;
            }
            dataLineOffset += dataLineStride;
            maskLineOffset += maskLineStride;
        }
    }

    public void accumulateDataFloat(PixelAccessor dataAccessor, Raster dataTile, PixelAccessor maskAccessor,
                                    Raster maskTile, Rectangle r, String unit) {
        final int[] bins = this.bins;
        final double lowValue = this.lowValue;
        final double highValue = this.highValue;
        final double binWidth = this.binWidth;

        final UnpackedImageData duid = dataAccessor.getPixels(dataTile, r, DataBuffer.TYPE_FLOAT, false);
        final float[] data = duid.getFloatData(0);
        final int dataPixelStride = duid.pixelStride;
        final int dataLineStride = duid.lineStride;
        final int dataBandOffset = duid.bandOffsets[0];

        byte[] mask = null;
        int maskPixelStride = 0;
        int maskLineStride = 0;
        int maskBandOffset = 0;
        if (maskAccessor != null) {
            UnpackedImageData muid = maskAccessor.getPixels(maskTile, r, DataBuffer.TYPE_BYTE, false);
            mask = muid.getByteData(0);
            maskPixelStride = muid.pixelStride;
            maskLineStride = muid.lineStride;
            maskBandOffset = muid.bandOffsets[0];
        }

        final int width = r.width;
        final int height = r.height;

        int dataLineOffset = dataBandOffset;
        int maskLineOffset = maskBandOffset;
        for (int y = 0; y < height; y++) {
            int dataPixelOffset = dataLineOffset;
            int maskPixelOffset = maskLineOffset;
            for (int x = 0; x < width; x++) {
                if (mask == null || mask[maskPixelOffset] != 0) {
                    double d = data[dataPixelOffset];
                    if (d >= lowValue && d <= highValue) {
                        int i = (int) ((d - lowValue) / binWidth);
                        i = i == bins.length ? i - 1 : i;
                        bins[i]++;
                    }
                }
                dataPixelOffset += dataPixelStride;
                maskPixelOffset += maskPixelStride;
            }
            dataLineOffset += dataLineStride;
            maskLineOffset += maskLineStride;
        }
    }

    public void accumulateDataDouble(PixelAccessor dataAccessor, Raster dataTile, PixelAccessor maskAccessor,
                                     Raster maskTile, Rectangle r, String unit) {
        final int[] bins = this.bins;
        final double lowValue = this.lowValue;
        final double highValue = this.highValue;
        final double binWidth = this.binWidth;

        final UnpackedImageData duid = dataAccessor.getPixels(dataTile, r, DataBuffer.TYPE_DOUBLE, false);
        final double[] data = duid.getDoubleData(0);
        final int dataPixelStride = duid.pixelStride;
        final int dataLineStride = duid.lineStride;
        final int dataBandOffset = duid.bandOffsets[0];

        byte[] mask = null;
        int maskPixelStride = 0;
        int maskLineStride = 0;
        int maskBandOffset = 0;
        if (maskAccessor != null) {
            UnpackedImageData muid = maskAccessor.getPixels(maskTile, r, DataBuffer.TYPE_BYTE, false);
            mask = muid.getByteData(0);
            maskPixelStride = muid.pixelStride;
            maskLineStride = muid.lineStride;
            maskBandOffset = muid.bandOffsets[0];
        }

        final int width = r.width;
        final int height = r.height;

        int dataLineOffset = dataBandOffset;
        int maskLineOffset = maskBandOffset;
        for (int y = 0; y < height; y++) {
            int dataPixelOffset = dataLineOffset;
            int maskPixelOffset = maskLineOffset;
            for (int x = 0; x < width; x++) {
                if (mask == null || mask[maskPixelOffset] != 0) {
                    double d = data[dataPixelOffset];
                    if (d >= lowValue && d <= highValue) {
                        int i = (int) ((d - lowValue) / binWidth);
                        i = i == bins.length ? i - 1 : i;
                        bins[i]++;
                    }
                }
                dataPixelOffset += dataPixelStride;
                maskPixelOffset += maskPixelStride;
            }
            dataLineOffset += dataLineStride;
            maskLineOffset += maskLineStride;
        }
    }
}
