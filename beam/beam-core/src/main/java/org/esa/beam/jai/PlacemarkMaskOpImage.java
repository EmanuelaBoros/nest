package org.esa.beam.jai;

import org.esa.beam.framework.datamodel.*;

import javax.media.jai.PlanarImage;
import javax.media.jai.RasterFactory;
import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.image.*;

/**
 * Creates a mask image for a given {@link org.esa.beam.framework.datamodel.RasterDataNode}.
 * The resulting image will have a single-band, interleaved sample model
 * with sample values 255 or 0.
 */
public class PlacemarkMaskOpImage extends SingleBandedOpImage {
    private static final byte FALSE = (byte) 0;
    private static final byte TRUE = (byte) 255;
    private final ColorModel colorModel;
    private final Product product;
    private final PlacemarkDescriptor placemarkDescriptor;
    private final int placemarkSize;

    public PlacemarkMaskOpImage(Product product,
                                PlacemarkDescriptor placemarkDescriptor,
                                int placemarkSize,
                                int width,
                                int height,
                                ResolutionLevel level) {
        super(DataBuffer.TYPE_BYTE,
              width,
              height,
              product.getPreferredTileSize(),
              null,
              level);
        this.product = product;
        this.placemarkDescriptor = placemarkDescriptor;
        this.placemarkSize = placemarkSize;
        this.colorModel = createColorModel(getSampleModel());
    }

    @Override
    protected void computeRect(PlanarImage[] sourceImages, WritableRaster tile, Rectangle destRect) {
        final BufferedImage image = new BufferedImage(colorModel, RasterFactory.createWritableRaster(tile.getSampleModel(), tile.getDataBuffer(), new Point(0, 0)), false, null);
        final Graphics2D graphics2D = image.createGraphics();
        graphics2D.translate(-tile.getMinX(), -tile.getMinY());
        graphics2D.setColor(Color.WHITE);

        ProductNodeGroup<Pin> pinGroup = getPlacemarkGroup();
        Pin[] placemarks = pinGroup.toArray(new Pin[pinGroup.getNodeCount()]);
        for (Pin placemark : placemarks) {
            final PixelPos pixelPos = placemark.getPixelPos();
            if (pixelPos != null) {
                final int x = (int) pixelPos.x - placemarkSize / 2;
                final int y = (int) pixelPos.y - placemarkSize / 2;
                graphics2D.fillRect(x, y, placemarkSize, placemarkSize);
            }
        }
        graphics2D.dispose();

        final byte[] data = ((DataBufferByte) tile.getDataBuffer()).getData();
        for (int i = 0; i < data.length; i++) {
            data[i] = (data[i] != 0) ? TRUE : FALSE;
        }
    }

    private ProductNodeGroup<Pin> getPlacemarkGroup() {
        return placemarkDescriptor.getPlacemarkGroup(product);
    }
}