package org.esa.beam.dataio.dimap.spi;

import com.bc.ceres.binding.PropertyContainer;
import static org.esa.beam.dataio.dimap.DimapProductConstants.*;
import org.esa.beam.framework.datamodel.Mask;
import org.esa.beam.framework.datamodel.Product;
import org.jdom.Document;
import org.jdom.Element;
import org.jdom.JDOMException;
import org.jdom.input.SAXBuilder;
import static org.junit.Assert.*;
import org.junit.Test;

import java.awt.Color;
import java.io.IOException;
import java.io.InputStream;

public class RangeTypePersistableTest {


    /*
    <Mask type="Range">
      <NAME value="myRange" />
      <DESCRIPTION value="Carefully defined range" />
      <COLOR red="0" green="255" blue="0" alpha="128" />
      <TRANSPARENCY value="0.78" />
      <MINIMUM value="0.35" />
      <MAXIMUM value="0.76" />
      <RASTER value="reflectance_13" />
    </Mask>
    */
    @Test
    public void createXmlFromObject() {
        final Mask.RangeType rangeType = new Mask.RangeType();
        final Mask mask = new Mask("myRange", 10, 10, rangeType);
        mask.setDescription("Carefully defined range");
        mask.setImageColor(new Color(0, 255, 0, 128));
        mask.setImageTransparency(0.78);
        final PropertyContainer config = mask.getImageConfig();
        config.setValue(Mask.RangeType.PROPERTY_NAME_MINIMUM, 0.35);
        config.setValue(Mask.RangeType.PROPERTY_NAME_MAXIMUM, 0.76);
        config.setValue(Mask.RangeType.PROPERTY_NAME_RASTER, "reflectance_13");

        final RangeTypePersistable persistable = new RangeTypePersistable();
        final Element element = persistable.createXmlFromObject(mask);
        assertNotNull(element);
        assertEquals(TAG_MASK, element.getName());
        assertEquals(Mask.RangeType.TYPE_NAME, getAttributeString(element, ATTRIB_TYPE));

        final Element name = element.getChild(TAG_NAME);
        assertEquals("myRange", getAttributeString(name, ATTRIB_VALUE));

        final Element description = element.getChild(TAG_DESCRIPTION);
        assertEquals("Carefully defined range", getAttributeString(description, ATTRIB_VALUE));

        final Element color = element.getChild(TAG_COLOR);
        assertEquals(0, getAttributeInt(color, ATTRIB_RED));
        assertEquals(255, getAttributeInt(color, ATTRIB_GREEN));
        assertEquals(0, getAttributeInt(color, ATTRIB_BLUE));
        assertEquals(128, getAttributeInt(color, ATTRIB_ALPHA));

        final Element transparency = element.getChild(TAG_TRANSPARENCY);
        assertEquals(0.78, getAttributeDouble(transparency, ATTRIB_VALUE), 0.0);


        final Element minimum = element.getChild(TAG_MINIMUM);
        assertEquals(0.35, getAttributeDouble(minimum, ATTRIB_VALUE), 0.0);
        final Element maximum = element.getChild(TAG_MAXIMUM);
        assertEquals(0.76, getAttributeDouble(maximum, ATTRIB_VALUE), 0.0);
        final Element raster = element.getChild(TAG_RASTER);
        assertEquals("reflectance_13", getAttributeString(raster, ATTRIB_VALUE));
    }

    @Test
    public void createMaskFromXml() throws IOException, JDOMException {
        final DimapPersistable persistable = new RangeTypePersistable();
        final InputStream resourceStream = getClass().getResourceAsStream("RangeMask.xml");
        final Document document = new SAXBuilder().build(resourceStream);
        final Product product = new Product("P", "T", 10, 10);
        final Mask maskFromXml = (Mask) persistable.createObjectFromXml(document.getRootElement(), product);

        assertNotNull(maskFromXml);
        assertEquals(Mask.RangeType.class, maskFromXml.getImageType().getClass());
        assertEquals("myRange", maskFromXml.getName());
        assertEquals("Carefully defined range", maskFromXml.getDescription());
        assertEquals(0.78, maskFromXml.getImageTransparency(), 0.0);
        assertEquals(new Color(0, 255, 0, 128), maskFromXml.getImageColor());

        assertEquals(0.35, maskFromXml.getImageConfig().getValue(Mask.RangeType.PROPERTY_NAME_MINIMUM));
        assertEquals(0.76, maskFromXml.getImageConfig().getValue(Mask.RangeType.PROPERTY_NAME_MAXIMUM));
        assertEquals("reflectance_13", maskFromXml.getImageConfig().getValue(Mask.RangeType.PROPERTY_NAME_RASTER));
    }

    private int getAttributeInt(Element element, String attribName) {
        return Integer.parseInt(element.getAttribute(attribName).getValue());
    }

    private double getAttributeDouble(Element element, String attribName) {
        return Double.parseDouble(element.getAttribute(attribName).getValue());
    }

    private String getAttributeString(Element element, String attribName) {
        return element.getAttribute(attribName).getValue();
    }
}
