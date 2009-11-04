package org.esa.beam.jai;

import org.esa.beam.framework.datamodel.VectorData;
import org.opengis.feature.simple.SimpleFeatureType;
import org.opengis.feature.simple.SimpleFeature;
import org.opengis.referencing.operation.TransformException;
import org.opengis.referencing.operation.MathTransform2D;
import org.opengis.referencing.FactoryException;
import org.geotools.feature.FeatureCollection;
import org.geotools.feature.FeatureIterator;
import org.geotools.geometry.jts.LiteShape2;
import org.geotools.geometry.jts.Decimator;
import org.geotools.referencing.operation.transform.AffineTransform2D;

import javax.media.jai.PlanarImage;
import javax.media.jai.RasterFactory;
import java.awt.image.DataBuffer;
import java.awt.image.WritableRaster;
import java.awt.image.BufferedImage;
import java.awt.image.DataBufferByte;
import java.awt.Rectangle;
import java.awt.Point;
import java.awt.Graphics2D;
import java.awt.Color;
import java.awt.geom.AffineTransform;
import java.awt.geom.Point2D;

import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.Lineal;
import com.vividsolutions.jts.geom.Polygonal;
import com.vividsolutions.jts.geom.GeometryCollection;
import com.vividsolutions.jts.geom.Puntal;
import com.vividsolutions.jts.geom.Coordinate;

public class VectorDataMaskOpImage extends SingleBandedOpImage {
    private static final byte FALSE = (byte) 0;
    private static final byte TRUE = (byte) 255;
    private final VectorData vectorData;

    public VectorDataMaskOpImage(VectorData vectorData, ResolutionLevel level) {
        super(DataBuffer.TYPE_BYTE,
              vectorData.getProduct().getSceneRasterWidth(),
              vectorData.getProduct().getSceneRasterHeight(),
              vectorData.getProduct().getPreferredTileSize(),
              null,
              level);
        this.vectorData = vectorData;
    }

    public VectorData getVectorData() {
        return vectorData;
    }

    @Override
    protected void computeRect(PlanarImage[] sourceImages, WritableRaster tile, Rectangle destRect) {
        final BufferedImage image = new BufferedImage(colorModel, RasterFactory.createWritableRaster(tile.getSampleModel(), tile.getDataBuffer(), new Point(0, 0)), false, null);
        final Graphics2D graphics2D = image.createGraphics();
        graphics2D.translate(-tile.getMinX(), -tile.getMinY());
        graphics2D.setColor(Color.WHITE);

        FeatureCollection<SimpleFeatureType,SimpleFeature> features = vectorData.getFeatureCollection();
        FeatureIterator<SimpleFeature> featureIterator = features.features();
        AffineTransform2D transform = new AffineTransform2D(AffineTransform.getScaleInstance(1.0 / getScale(), 1.0 / getScale()));
        while (featureIterator.hasNext()) {
            SimpleFeature feature = featureIterator.next();
            Object value = feature.getDefaultGeometry();
            if (value instanceof Geometry) {
                try {
                    renderGeometry((Geometry) value, graphics2D, transform);
                } catch (Exception ignored) {
                    // ignore
                }
            }
        }

        graphics2D.dispose();

        final byte[] data = ((DataBufferByte) tile.getDataBuffer()).getData();
        for (int i = 0; i < data.length; i++) {
            data[i] = (data[i] != 0) ? TRUE : FALSE;
        }
    }

    private static void renderGeometry(Geometry geom, Graphics2D graphics, MathTransform2D transform) throws Exception {
        if (geom instanceof Puntal) {
            Coordinate c = geom.getCoordinate();
            Point2D.Double pt = new Point2D.Double(c.x, c.y);
            transform.transform(pt, pt);
            graphics.drawLine((int)pt.x, (int)pt.y, (int)pt.x, (int)pt.y);
        } else if (geom instanceof Lineal) {
            LiteShape2 shape = new LiteShape2(geom, transform, null, false, true);
            graphics.draw(shape);
        } else if (geom instanceof Polygonal) {
            LiteShape2 shape = new LiteShape2(geom, transform, null, false, true);
            graphics.fill(shape);
        } else if (geom instanceof GeometryCollection) {
            GeometryCollection collection = (GeometryCollection) geom;
            for (int i = 0; i < collection.getNumGeometries(); i++) {
                renderGeometry(collection.getGeometryN(i), graphics, transform);
            }
        }
    }
}
