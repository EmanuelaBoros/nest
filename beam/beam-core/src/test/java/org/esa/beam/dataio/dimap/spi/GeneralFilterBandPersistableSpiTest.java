package org.esa.beam.dataio.dimap.spi;
/*
 * $Id: GeneralFilterBandPersistableSpiTest.java,v 1.2 2010-03-31 13:59:55 lveci Exp $
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

import junit.framework.TestCase;
import org.esa.beam.dataio.dimap.DimapProductConstants;
import org.esa.beam.framework.datamodel.Band;
import org.esa.beam.framework.datamodel.GeneralFilterBand;
import org.esa.beam.framework.datamodel.ProductData;
import org.jdom.Attribute;
import org.jdom.Element;

import java.util.ArrayList;

/**
 * Created by Marco Peters.
 *
 * @author Marco Peters
 * @version $Revision: 1.2 $ $Date: 2010-03-31 13:59:55 $
 */
public class GeneralFilterBandPersistableSpiTest extends TestCase {

    private GeneralFilterBandPersistableSpi _persistableSpi;

    @Override
    public void setUp() {
        _persistableSpi = new GeneralFilterBandPersistableSpi();
    }

    @Override
    protected void tearDown() throws Exception {
        _persistableSpi = null;
    }

    public void testCanDecode_GoodElement() {

        final Element bandInfo = new Element(DimapProductConstants.TAG_SPECTRAL_BAND_INFO);
        final Element filterInfo = new Element(DimapProductConstants.TAG_FILTER_BAND_INFO);
        final Attribute bandType = new Attribute(DimapProductConstants.ATTRIB_BAND_TYPE, "GeneralFilterBand");
        filterInfo.setAttribute(bandType);
        bandInfo.setContent(filterInfo);

        assertTrue(_persistableSpi.canDecode(bandInfo));
    }

    public void testCanDecode_NotSpectralBandInfo() {

        final Element element = new Element("SomeWhat");

        assertFalse(_persistableSpi.canDecode(element));
    }

    public void testCanDecode_NoBandType() {

        final Element bandInfo = new Element(DimapProductConstants.TAG_SPECTRAL_BAND_INFO);
        final Element filterInfo = new Element(DimapProductConstants.TAG_FILTER_BAND_INFO);
        bandInfo.setContent(filterInfo);

        assertFalse(_persistableSpi.canDecode(bandInfo));
    }

    public void testCanDecode_NotCorrectBandType() {
        final Element bandInfo = new Element(DimapProductConstants.TAG_SPECTRAL_BAND_INFO);
        final Element filterInfo = new Element(DimapProductConstants.TAG_FILTER_BAND_INFO);
        final Attribute bandType = new Attribute(DimapProductConstants.ATTRIB_BAND_TYPE, "VirtualBand");
        filterInfo.setAttribute(bandType);
        bandInfo.setContent(filterInfo);

        assertFalse(_persistableSpi.canDecode(bandInfo));
    }


    public void testCanPersist() {
        final Band source = new Band("b", ProductData.TYPE_INT8, 2, 2);
        final GeneralFilterBand gfb = new GeneralFilterBand("test", source, 3, 3, GeneralFilterBand.MAX);

        assertTrue(_persistableSpi.canPersist(gfb));

        assertFalse(_persistableSpi.canPersist(new ArrayList()));
        assertFalse(_persistableSpi.canPersist(new Object()));
        assertFalse(_persistableSpi.canPersist(new Band("b", ProductData.TYPE_INT8, 2, 2)));
    }

    public void testCreatePersistable() {
        assertNotNull(_persistableSpi.createPersistable());
    }
}
