/*
 * $Id: Pin.java,v 1.7 2010-02-08 21:57:50 lveci Exp $
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

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.Point;
import org.esa.beam.dataio.dimap.DimapProductConstants;
import org.esa.beam.dataio.dimap.DimapProductHelpers;
import org.esa.beam.framework.dataio.ProductSubsetDef;
import org.esa.beam.jai.ImageManager;
import org.esa.beam.util.Debug;
import org.esa.beam.util.Guardian;
import org.esa.beam.util.ObjectUtils;
import org.esa.beam.util.XmlWriter;
import org.geotools.feature.simple.SimpleFeatureTypeBuilder;
import org.jdom.Element;
import org.opengis.feature.simple.SimpleFeature;
import org.opengis.feature.simple.SimpleFeatureType;
import org.opengis.feature.type.AttributeDescriptor;

import java.awt.Color;
import java.awt.Rectangle;
import java.awt.geom.AffineTransform;
import java.awt.geom.NoninvertibleTransformException;
import java.awt.geom.Point2D;
import java.text.MessageFormat;
import java.util.List;

// todo - rename to Placemark (se - 20090126) 

/**
 * This class represents a pin.
 * <p/>
 * Pins are displayed as symbols at the image's pixel position corresponding to their geographical position. The name is
 * displayed as label next to the symbol. If the user moves the mouse over a pin, the textual description property shall
 * appear as tool-tip text. Single pins can be selected either by mouse-click or by the ? Prev./Next Pin tool. Pins are
 * contained in the active product and stored in DIMAP format. To share pins between products, the pins of a product can
 * be imported and exported.
 *
 * @author Sabine Embacher
 * @version $Revision: 1.7 $ $Date: 2010-02-08 21:57:50 $
 */
public class Pin extends ProductNode {

    public static final String PLACEMARK_FEATURE_TYPE_NAME = "Placemark";

    public static final String PROPERTY_NAME_LABEL = "label";
    public static final String PROPERTY_NAME_PIXELPOS = "pixelPos";
    public static final String PROPERTY_NAME_GEOPOS = "geoPos";
    public static final String PROPERTY_NAME_PINSYMBOL = "pinSymbol";

    private SimpleFeature feature;

    /**
     * Returns the type of features underlying all pins.
     *
     * @return the type of features underlying all pins.
     * @since BEAM 4.7
     */
    public static SimpleFeatureType getFeatureType() {
        return Holder.PLACEMARK_FEATURE_TYPE;
    }

    /**
     * Creates a new pin.
     *
     * @param name     the pin's name.
     * @param label    the pin's label.
     * @param pixelPos the pin's pixel position.
     * @deprecated since 4.1, use {@link Pin#Pin(String, String, String, PixelPos, GeoPos, PlacemarkSymbol)}
     */
    @Deprecated
    public Pin(String name, String label, PixelPos pixelPos) {
        this(name, label, "", pixelPos, null, PlacemarkSymbol.createDefaultPinSymbol());
    }

    /**
     * Creates a new pin.
     *
     * @param name   the pin's name.
     * @param label  the pin's label.
     * @param geoPos the pin's geo-position.
     * @deprecated since 4.1, use {@link Pin#Pin(String, String, String, PixelPos, GeoPos, PlacemarkSymbol)}
     */
    @Deprecated
    public Pin(String name, String label, GeoPos geoPos) {
        this(name, label, "", null, geoPos, PlacemarkSymbol.createDefaultPinSymbol());
    }

    /**
     * Creates a new pin.
     *
     * @param name        the pin's name.
     * @param label       the pin's label.
     * @param description the pin's description
     * @param pixelPos    the pin's pixel position
     * @param geoPos      the pin's geo-position.
     * @param symbol      the pin's symbol.
     * @deprecated since 4.7, use {@link Pin#Pin(String, String, String, PixelPos, GeoPos, PlacemarkSymbol, GeoCoding)}
     */
    @Deprecated
    public Pin(String name, String label, String description, PixelPos pixelPos, GeoPos geoPos,
               PlacemarkSymbol symbol) {
        this(name, label, description, pixelPos, geoPos, symbol, null);
    }

    /**
     * Creates a new pin.
     *
     * @param name        the pin's name.
     * @param label       the pin's label.
     * @param description the pin's description
     * @param pixelPos    the pin's pixel position
     * @param geoPos      the pin's geo-position.
     * @param symbol      the pin's symbol.
     * @param geoCoding   the pin's geo-coding.
     */
    public Pin(String name, String label, String description, PixelPos pixelPos, GeoPos geoPos,
               PlacemarkSymbol symbol, GeoCoding geoCoding) {
        super(name, description);
        if (pixelPos == null && geoPos == null) {
            throw new IllegalArgumentException("pixelPos == null && geoPos == null");
        }
        // todo: rq/?? - we need to ensure that the pixel position can be calculated, though many tests fail then
        // if (pixelPos == null && geoCoding == null || !geoCoding.canGetPixelPos()) {
        //     throw new IllegalArgumentException("pixelPos == null && geoCoding == null || !geoCoding.canGetPixelPos()");
        // }

        feature = createFeature(name, label, pixelPos, geoPos, symbol, geoCoding);
    }

    Pin(SimpleFeature feature) {
        super(feature.getID(), "");
        this.feature = feature;
    }

    /**
     * Returns the {@link SimpleFeature}, underlying this pin.
     *
     * @return the {@link SimpleFeature} underlying this pin.
     * @since BEAM 4.7
     */
    public final SimpleFeature getFeature() {
        return feature;
    }

    /**
     * Sets this pin's label.
     *
     * @param label the label, if {@code null} an empty label is set.
     */
    public void setLabel(String label) {
        if (label == null) {
            label = "";
        }
        if (!label.equals(feature.getAttribute(PROPERTY_NAME_LABEL))) {
            feature.setAttribute(PROPERTY_NAME_LABEL, label);
            fireProductNodeChanged(PROPERTY_NAME_LABEL);
        }
    }

    /**
     * Returns this pin's label.
     *
     * @return the label, cannot be {@code null}.
     */
    public String getLabel() {
        return (String) feature.getAttribute(PROPERTY_NAME_LABEL);
    }

    /**
     * Returns an estimated, raw storage size in bytes of this pin.
     *
     * @param subsetDef if not {@code null} the subset may limit the size returned.
     * @return the estimated size in bytes.
     */
    @Override
    public long getRawStorageSize(ProductSubsetDef subsetDef) {
        return 256;
    }

    /**
     * Accepts the given visitor. This method implements the well known 'Visitor' design pattern of the gang-of-four.
     * The visitor pattern allows to define new operations on the product data model without the need to add more code
     * to it. The new operation is implemented by the visitor.
     *
     * @param visitor the visitor
     */
    @Override
    public void acceptVisitor(ProductVisitor visitor) {
    }

    public PlacemarkSymbol getSymbol() {
        return (PlacemarkSymbol) feature.getAttribute("symbol");
    }

    public void setSymbol(final PlacemarkSymbol symbol) {
        Guardian.assertNotNull("symbol", symbol);
        if (getSymbol() != symbol) {
            feature.setAttribute("symbol", symbol);
            fireProductNodeChanged(PROPERTY_NAME_PINSYMBOL);
        }
    }

    public void setPixelPos(PixelPos pixelPos) {
        setPixelCoordinate(pixelPos, true);
    }

    /**
     * Returns this pin's pixel position.
     *
     * @return this pin's pixel position. If this pin's pixel position is {@code null}, the pixel
     *         position returned is calculated from this pin's geo-position, if possible.
     */
    public PixelPos getPixelPos() {
        PixelPos pixelPos = toPixelPos(getPixelCoordinate());
        if (pixelPos == null && canComputePixelPos()) {
            final GeoPos geoPos = toGeoPos(getGeoCoordinate());
            if (geoPos != null) {
                pixelPos = getProduct().getGeoCoding().getPixelPos(geoPos, null);
            }
        }
        if (pixelPos == null) {
            return null;
        }
        if (getProduct() != null) {
            final int w = getProduct().getSceneRasterWidth();
            final int h = getProduct().getSceneRasterHeight();
            final Rectangle bounds = new Rectangle(0, 0, w, h);
            if (!bounds.contains(pixelPos.x, pixelPos.y)) {
                return null;
            }
        }

        return pixelPos;
    }

    public void setGeoPos(GeoPos geoPos) {
        setGeoCoordinate(geoPos);
    }

    public GeoPos getGeoPos() {
        GeoPos geoPos = toGeoPos(getGeoCoordinate());
        if (geoPos == null && canComputeGeoPos()) {
            final PixelPos pixelPos = toPixelPos(getPixelCoordinate());
            if (pixelPos != null) {
                geoPos = getProduct().getGeoCoding().getGeoPos(pixelPos, null);
            }
        }
        if (geoPos == null) {
            return null;
        }
        return geoPos;
    }

    public void writeXML(XmlWriter writer, int indent) {
        Guardian.assertNotNull("writer", writer);
        Guardian.assertGreaterThan("indent", indent, -1);

        final String[][] attributes = {new String[]{DimapProductConstants.ATTRIB_NAME, getName()}};
        final String[] pinTags = XmlWriter.createTags(indent, DimapProductConstants.TAG_PLACEMARK, attributes);
        writer.println(pinTags[0]);
        indent++;
        writer.printLine(indent, DimapProductConstants.TAG_PLACEMARK_LABEL, getLabel());
        writer.printLine(indent, DimapProductConstants.TAG_PLACEMARK_DESCRIPTION, getDescription());
        final GeoPos geoPos = getGeoPos();
        if (geoPos != null) {
            writer.printLine(indent, DimapProductConstants.TAG_PLACEMARK_LATITUDE, geoPos.lat);
            writer.printLine(indent, DimapProductConstants.TAG_PLACEMARK_LONGITUDE, geoPos.lon);
        }
        final PixelPos pixelPos = getPixelPos();
        if (pixelPos != null) {
            writer.printLine(indent, DimapProductConstants.TAG_PLACEMARK_PIXEL_X, pixelPos.x);
            writer.printLine(indent, DimapProductConstants.TAG_PLACEMARK_PIXEL_Y, pixelPos.y);
        }
        final Color fillColor = (Color) getSymbol().getFillPaint();
        if (fillColor != null) {
            writeColor(DimapProductConstants.TAG_PLACEMARK_FILL_COLOR, indent, fillColor, writer);
        }
        final Color outlineColor = getSymbol().getOutlineColor();
        if (outlineColor != null) {
            writeColor(DimapProductConstants.TAG_PLACEMARK_OUTLINE_COLOR, indent, outlineColor, writer);
        }
        writer.println(pinTags[1]);
    }

    // todo - move this method into a new DimapPersistable

    private void writeColor(final String tagName, final int indent, final Color color, final XmlWriter writer) {
        final String[] colorTags = XmlWriter.createTags(indent, tagName);
        writer.println(colorTags[0]);
        DimapProductHelpers.printColorTag(indent + 1, color, writer);
        writer.println(colorTags[1]);
    }

    /**
     * Creates a new GCP from an XML element.
     *
     * @param element the element.
     * @return the GCP created.
     * @throws NullPointerException     if element is null
     * @throws IllegalArgumentException if element is invalid
     * @deprecated since BEAM 4.7, use {@link #createPlacemark(Element, PlacemarkSymbol, GeoCoding)} instead
     */
    @Deprecated
    public static Pin createGcp(Element element) {
        return createPlacemark(element, PlacemarkSymbol.createDefaultGcpSymbol(), null);
    }

    /**
     * Creates a new pin from an XML element.
     *
     * @param element the element.
     * @return the pin created.
     * @throws NullPointerException     if element is null
     * @throws IllegalArgumentException if element is invalid
     * @deprecated since BEAM 4.7, use {@link #createPlacemark(org.jdom.Element, PlacemarkSymbol, GeoCoding)} instead
     */
    @Deprecated
    public static Pin createPin(Element element) {
        return createPlacemark(element, PlacemarkSymbol.createDefaultPinSymbol(), null);
    }

    /**
     * Creates a new pin from an XML element and a given symbol.
     *
     * @param element the element.
     * @param symbol  the symbol.
     * @return the pin created.
     * @throws NullPointerException     if element is null
     * @throws IllegalArgumentException if element is invalid
     * @deprecated since BEAM 4.7, use {@link #createPlacemark(Element, PlacemarkSymbol, GeoCoding)} instead
     */
    @Deprecated
    public static Pin createPlacemark(Element element, PlacemarkSymbol symbol) {
        return createPlacemark(element, symbol, null);
    }
    /**
     * Creates a new pin from an XML element and a given symbol.
     *
     * @param element the element.
     * @param symbol  the symbol.
     * @param geoCoding the geoCoding to used by the placemark. Can be <code>null</code>.
     *
     * @return the pin created.
     *
     * @throws NullPointerException     if element is null
     * @throws IllegalArgumentException if element is invalid
     */
    public static Pin createPlacemark(Element element, PlacemarkSymbol symbol, GeoCoding geoCoding) {
        if (!DimapProductConstants.TAG_PLACEMARK.equals(element.getName()) &&
                !DimapProductConstants.TAG_PIN.equals(element.getName())) {
            throw new IllegalArgumentException(MessageFormat.format("Element ''{0}'' or ''{1}'' expected.",
                                                                    DimapProductConstants.TAG_PLACEMARK,
                                                                    DimapProductConstants.TAG_PIN));
        }
        final String name = element.getAttributeValue(DimapProductConstants.ATTRIB_NAME);
        if (name == null) {
            throw new IllegalArgumentException(MessageFormat.format("Missing attribute ''{0}''.",
                                                                    DimapProductConstants.ATTRIB_NAME));
        }
        if (!Pin.isValidNodeName(name)) {
            throw new IllegalArgumentException(MessageFormat.format("Invalid placemark name ''{0}''.", name));
        }

        String label = element.getChildTextTrim(DimapProductConstants.TAG_PLACEMARK_LABEL);
        if (label == null) {
            label = name;
        }
        final String description = element.getChildTextTrim(DimapProductConstants.TAG_PLACEMARK_DESCRIPTION);
        final String latText = element.getChildTextTrim(DimapProductConstants.TAG_PLACEMARK_LATITUDE);
        final String lonText = element.getChildTextTrim(DimapProductConstants.TAG_PLACEMARK_LONGITUDE);
        final String posXText = element.getChildTextTrim(DimapProductConstants.TAG_PLACEMARK_PIXEL_X);
        final String posYText = element.getChildTextTrim(DimapProductConstants.TAG_PLACEMARK_PIXEL_Y);

        GeoPos geoPos = null;
        if (latText != null && lonText != null) {
            try {
                float lat = Float.parseFloat(latText);
                float lon = Float.parseFloat(lonText);
                geoPos = new GeoPos(lat, lon);
            } catch (NumberFormatException e) {
                throw new IllegalArgumentException("Invalid geo-position.");
            }
        }
        PixelPos pixelPos = null;
        if (posXText != null && posYText != null) {
            try {
                float pixelX = Float.parseFloat(posXText);
                float pixelY = Float.parseFloat(posYText);
                pixelPos = new PixelPos(pixelX, pixelY);
            } catch (NumberFormatException e) {
                throw new IllegalArgumentException("Invalid pixel-position.");
            }
        }
        if (geoPos == null && pixelPos == null) {
            throw new IllegalArgumentException("Neither geo-position nor pixel-position given.");
        }
        final Color fillColor = createColor(element.getChild(DimapProductConstants.TAG_PLACEMARK_FILL_COLOR));
        if (fillColor != null) {
            symbol.setFillPaint(fillColor);
        }
        final Color outlineColor = createColor(element.getChild(DimapProductConstants.TAG_PLACEMARK_OUTLINE_COLOR));
        if (outlineColor != null) {
            symbol.setOutlineColor(outlineColor);
        }

        return new Pin(name, label, description, pixelPos, geoPos, symbol, geoCoding);
    }

    // todo - move this method into a new DimapPersistable

    private static Color createColor(Element elem) {
        if (elem != null) {
            Element colorElem = elem.getChild(DimapProductConstants.TAG_COLOR);
            if (colorElem != null) {
                try {
                    return DimapProductHelpers.createColor(colorElem);
                } catch (NumberFormatException e) {
                    Debug.trace(e);
                } catch (IllegalArgumentException e) {
                    Debug.trace(e);
                }
            }
        }
        return null;
    }

    /**
     * Releases all of the resources used by this object instance and all of its owned children. Its primary use is to
     * allow the garbage collector to perform a vanilla job.
     * <p/>
     * <p>This method should be called only if it is for sure that this object instance will never be used again. The
     * results of referencing an instance of this class after a call to <code>dispose()</code> are undefined.
     * <p/>
     * <p>Overrides of this method should always call <code>super.dispose();</code> after disposing this instance.
     */
    @Override
    public void dispose() {
        if (feature != null) {
            final PlacemarkSymbol symbol = getSymbol();
            if (symbol != null) {
                symbol.dispose();
            }
            for (final AttributeDescriptor attributeDescriptor : feature.getFeatureType().getAttributeDescriptors()) {
                feature.setAttribute(attributeDescriptor.getLocalName(), null);
            }
            feature = null;
        }
        super.dispose();
    }


    private boolean canComputeGeoPos() {
        return (getProduct() != null
                && getProduct().getGeoCoding() != null
                && getProduct().getGeoCoding().canGetGeoPos());
    }

    private boolean canComputePixelPos() {
        return (getProduct() != null
                && getProduct().getGeoCoding() != null
                && getProduct().getGeoCoding().canGetPixelPos());
    }

    public void updatePixelPos(GeoCoding geoCoding) {
        if (getGeoPos() != null && geoCoding != null && geoCoding.canGetPixelPos()) {
            setPixelPos(geoCoding.getPixelPos(getGeoPos(), null));
        }
    }

    public void updatePixelPos() {
        final Point point = (Point) feature.getDefaultGeometry();
        if (point != null) {
            if (getProduct() != null) {
                final AffineTransform i2m = ImageManager.getImageToModelTransform(getProduct().getGeoCoding());
                final PixelPos pixelPos = new PixelPos((float) point.getX(), (float) point.getY());
                try {
                    i2m.inverseTransform(pixelPos, pixelPos);
                } catch (NoninvertibleTransformException e) {
                    // ignore
                }
                setPixelCoordinate(pixelPos, false);
            }
        }
    }

    private void updateGeometry(PixelPos pixelPos) {
        if (getProduct() != null) {
            final AffineTransform i2m = ImageManager.getImageToModelTransform(getProduct().getGeoCoding());
            i2m.transform(pixelPos, pixelPos);

            final Point point = (Point) feature.getDefaultGeometry();
            point.getCoordinate().setCoordinate(toCoordinate(pixelPos));
            point.geometryChanged();
        }
    }

    private Coordinate getPixelCoordinate() {
        final Point point = (Point) feature.getAttribute(PROPERTY_NAME_PIXELPOS);
        if (point != null) {
            return point.getCoordinate();
        }
        return null;
    }

    private void setPixelCoordinate(PixelPos pixelPos, boolean updateGeometry) {
        final Coordinate newCoordinate = toCoordinate(pixelPos);
        final Coordinate oldCoordinate = getPixelCoordinate();
        if (!ObjectUtils.equalObjects(oldCoordinate, newCoordinate)) {
            if (oldCoordinate == null) {
                final GeometryFactory geometryFactory = new GeometryFactory();
                feature.setAttribute(PROPERTY_NAME_PIXELPOS, geometryFactory.createPoint(newCoordinate));
            } else {
                final Point point = (Point) feature.getAttribute(PROPERTY_NAME_PIXELPOS);
                point.getCoordinate().setCoordinate(newCoordinate);
                point.geometryChanged();
            }
            if (updateGeometry) {
                updateGeometry(pixelPos);
            }
            fireProductNodeChanged(PROPERTY_NAME_PIXELPOS);
        }
    }

    private Coordinate getGeoCoordinate() {
        final Point point = (Point) feature.getAttribute(PROPERTY_NAME_GEOPOS);
        if (point != null) {
            return point.getCoordinate();
        }
        return null;
    }

    private void setGeoCoordinate(GeoPos geoPos) {
        final Coordinate newCoordinate = toCoordinate(geoPos);
        final Coordinate oldCoordinate = getGeoCoordinate();
        if (!ObjectUtils.equalObjects(oldCoordinate, newCoordinate)) {
            if (oldCoordinate == null) {
                final GeometryFactory geometryFactory = new GeometryFactory();
                feature.setAttribute(PROPERTY_NAME_GEOPOS, geometryFactory.createPoint(newCoordinate));
            } else {
                final Point point = (Point) feature.getAttribute(PROPERTY_NAME_GEOPOS);
                point.getCoordinate().setCoordinate(newCoordinate);
                point.geometryChanged();
            }
            fireProductNodeChanged(PROPERTY_NAME_GEOPOS);
        }
    }

    private static SimpleFeatureType createFeatureType(String name) {
        final SimpleFeatureTypeBuilder builder = new SimpleFeatureTypeBuilder();

        final SimpleFeatureType superType = PlainFeatureFactory.createPlainFeatureType(name, Point.class, null);
        final List<AttributeDescriptor> list = superType.getAttributeDescriptors();
        for (AttributeDescriptor descriptor : list) {
            builder.add(descriptor);
        }
        builder.setName(name);
        builder.add(Pin.PROPERTY_NAME_LABEL, String.class);
        builder.add(Pin.PROPERTY_NAME_PIXELPOS, Point.class);
        builder.add(Pin.PROPERTY_NAME_GEOPOS, Point.class);
        builder.add("symbol", PlacemarkSymbol.class);

        return builder.buildFeatureType();
    }

    private static SimpleFeature createFeature(String name, String label, PixelPos pixelPos, GeoPos geoPos,
                                               PlacemarkSymbol symbol, GeoCoding geoCoding) {
        final SimpleFeatureType featureType = getFeatureType();
        final GeometryFactory geometryFactory = new GeometryFactory();
        final AffineTransform i2m = ImageManager.getImageToModelTransform(geoCoding);
        PixelPos imagePos = pixelPos;
        if (imagePos == null && geoCoding != null && geoCoding.canGetPixelPos()) {
            imagePos = geoCoding.getPixelPos(geoPos, imagePos);
        }
        if (imagePos == null) {
            imagePos = new PixelPos();
            imagePos.setInvalid();
        }
        final Point geometry = geometryFactory.createPoint(toCoordinate(i2m.transform(imagePos, null)));
        final SimpleFeature feature = PlainFeatureFactory.createPlainFeature(featureType, name, geometry, null);

        if (pixelPos != null) {
            feature.setAttribute(Pin.PROPERTY_NAME_PIXELPOS, geometryFactory.createPoint(toCoordinate(pixelPos)));
        }
        if (geoPos != null) {
            feature.setAttribute(Pin.PROPERTY_NAME_GEOPOS, geometryFactory.createPoint(toCoordinate(geoPos)));
        }
        if (label == null) {
            feature.setAttribute(Pin.PROPERTY_NAME_LABEL, "");
        } else {
            feature.setAttribute(Pin.PROPERTY_NAME_LABEL, label);
        }
        feature.setAttribute("symbol", symbol);

        return feature;
    }

    private static Coordinate toCoordinate(GeoPos geoPos) {
        if (geoPos != null) {
            return new Coordinate(geoPos.getLon(), geoPos.getLat());
        }
        return null;
    }

    private static Coordinate toCoordinate(Point2D pixelPos) {
        if (pixelPos != null) {
            return new Coordinate(pixelPos.getX(), pixelPos.getY());
        }
        return null;
    }

    private static GeoPos toGeoPos(Coordinate coordinate) {
        if (coordinate != null) {
            return new GeoPos((float) coordinate.y, (float) coordinate.x);
        }
        return null;
    }

    private static PixelPos toPixelPos(Coordinate coordinate) {
        if (coordinate != null) {
            return new PixelPos((float) coordinate.x, (float) coordinate.y);
        }
        return null;
    }

    private static class Holder {

        private static final SimpleFeatureType PLACEMARK_FEATURE_TYPE = createFeatureType(PLACEMARK_FEATURE_TYPE_NAME);
    }
}
