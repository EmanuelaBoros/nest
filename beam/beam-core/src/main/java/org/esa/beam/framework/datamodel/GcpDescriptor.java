package org.esa.beam.framework.datamodel;

import java.awt.Image;
import java.awt.Point;

/**
 * Created by Marco Peters.
 *
 * @author Marco Peters
 * @version $Revision: 1.3 $ $Date: 2009-12-23 16:42:11 $
 */
public class GcpDescriptor implements PlacemarkDescriptor {

    public final static GcpDescriptor INSTANCE = new GcpDescriptor();

    private GcpDescriptor() {
    }

    @Override
    public String getShowLayerCommandId() {
        return "showGcpOverlay";
    }

    @Override
    public String getRoleName() {
        return "gcp";
    }

    @Override
    public String getRoleLabel() {
        return "GCP";
    }

    @Override
    public Image getCursorImage() {
        return null;
    }

    @Override
    public Point getCursorHotSpot() {
        return new Point();
    }

    @Override
    public PlacemarkGroup getPlacemarkGroup(Product product) {
        return product.getGcpGroup();
    }

    @Override
    public PlacemarkSymbol createDefaultSymbol() {
        return PlacemarkSymbol.createDefaultGcpSymbol();
    }

    @Override
    public PixelPos updatePixelPos(GeoCoding geoCoding, GeoPos geoPos, PixelPos pixelPos) {
        if (geoCoding == null || !geoCoding.canGetPixelPos()) {
            return pixelPos;
        }
        return geoCoding.getPixelPos(geoPos, pixelPos);
    }

    @Override
    public GeoPos updateGeoPos(GeoCoding geoCoding, PixelPos pixelPos, GeoPos geoPos) {
        return geoPos;
    }
}
