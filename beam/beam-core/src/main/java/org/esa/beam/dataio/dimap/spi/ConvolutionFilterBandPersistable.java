/*
 * $Id: ConvolutionFilterBandPersistable.java,v 1.1 2009-04-28 14:39:32 lveci Exp $
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
package org.esa.beam.dataio.dimap.spi;

import org.esa.beam.dataio.dimap.DimapProductConstants;
import org.esa.beam.framework.datamodel.ConvolutionFilterBand;
import org.esa.beam.framework.datamodel.Kernel;
import org.esa.beam.framework.datamodel.Product;
import org.esa.beam.framework.datamodel.ProductData;
import org.esa.beam.framework.datamodel.RasterDataNode;
import org.esa.beam.util.StringUtils;
import org.jdom.Element;

import java.util.ArrayList;

/**
 * Created by Marco Peters.
 *
 * <p><i>Note that this class is not yet public API. Interface may chhange in future releases.</i></p>
 *
 * @author Marco Peters
 * @version $Revision: 1.1 $ $Date: 2009-04-28 14:39:32 $
 */
class ConvolutionFilterBandPersistable implements DimapPersistable {

    public Object createObjectFromXml(Element element, Product product) {
        final Element filterInfo = element.getChild(DimapProductConstants.TAG_FILTER_BAND_INFO);
        final Element kernelInfo = filterInfo.getChild(DimapProductConstants.TAG_FILTER_KERNEL);
        final String kernelDataString = kernelInfo.getChildTextTrim(DimapProductConstants.TAG_KERNEL_DATA);
        final double[] kernelData = StringUtils.toDoubleArray(kernelDataString, ",");
        final Kernel kernel = new Kernel(
                Integer.parseInt(kernelInfo.getChildTextTrim(DimapProductConstants.TAG_KERNEL_WIDTH)),
                Integer.parseInt(kernelInfo.getChildTextTrim(DimapProductConstants.TAG_KERNEL_HEIGHT)),
                Double.parseDouble(kernelInfo.getChildTextTrim(DimapProductConstants.TAG_KERNEL_FACTOR)),
                kernelData);
        final String sourceName = filterInfo.getChildTextTrim(DimapProductConstants.TAG_FILTER_SOURCE);
        final String bandName = element.getChildTextTrim(DimapProductConstants.TAG_BAND_NAME);
        final RasterDataNode sourceNode = product.getRasterDataNode(sourceName);
        final ConvolutionFilterBand cfb = new ConvolutionFilterBand(bandName, sourceNode,kernel);
        cfb.setDescription(element.getChildTextTrim(DimapProductConstants.TAG_BAND_DESCRIPTION));
        cfb.setUnit(element.getChildTextTrim(DimapProductConstants.TAG_PHYSICAL_UNIT));
        cfb.setSolarFlux(Float.parseFloat(element.getChildTextTrim(DimapProductConstants.TAG_SOLAR_FLUX)));
        cfb.setSpectralWavelength(Float.parseFloat(element.getChildTextTrim(DimapProductConstants.TAG_BAND_WAVELEN)));
        cfb.setSpectralBandwidth(Float.parseFloat(element.getChildTextTrim(DimapProductConstants.TAG_BANDWIDTH)));
        cfb.setScalingFactor(Double.parseDouble(element.getChildTextTrim(DimapProductConstants.TAG_SCALING_FACTOR)));
        cfb.setScalingOffset(Double.parseDouble(element.getChildTextTrim(DimapProductConstants.TAG_SCALING_OFFSET)));
        cfb.setLog10Scaled(Boolean.parseBoolean(element.getChildTextTrim(DimapProductConstants.TAG_SCALING_LOG_10)));
        cfb.setNoDataValueUsed(
                Boolean.parseBoolean(element.getChildTextTrim(DimapProductConstants.TAG_NO_DATA_VALUE_USED)));
        cfb.setNoDataValue(Double.parseDouble(element.getChildTextTrim(DimapProductConstants.TAG_NO_DATA_VALUE)));

        return cfb;
    }

    public Element createXmlFromObject(Object object) {
        final ConvolutionFilterBand cfb = (ConvolutionFilterBand) object;
        final ArrayList contentList = new ArrayList();
        contentList.add(createElement(DimapProductConstants.TAG_BAND_INDEX, String.valueOf(cfb.getProduct().getBandIndex(cfb.getName()))));
        contentList.add(createElement(DimapProductConstants.TAG_BAND_NAME, cfb.getName()));
        contentList.add(createElement(DimapProductConstants.TAG_BAND_DESCRIPTION, cfb.getDescription()));
        contentList.add(createElement(DimapProductConstants.TAG_DATA_TYPE, ProductData.getTypeString(cfb.getDataType())));
        contentList.add(createElement(DimapProductConstants.TAG_PHYSICAL_UNIT, cfb.getUnit()));
        contentList.add(createElement(DimapProductConstants.TAG_SOLAR_FLUX, String.valueOf(cfb.getSolarFlux())));
        contentList.add(createElement(DimapProductConstants.TAG_BAND_WAVELEN, String.valueOf(cfb.getSpectralWavelength())));
        contentList.add(createElement(DimapProductConstants.TAG_BANDWIDTH, String.valueOf(cfb.getSpectralBandwidth())));
        contentList.add(createElement(DimapProductConstants.TAG_SCALING_FACTOR, String.valueOf(cfb.getScalingFactor())));
        contentList.add(createElement(DimapProductConstants.TAG_SCALING_OFFSET, String.valueOf(cfb.getScalingOffset())));
        contentList.add(createElement(DimapProductConstants.TAG_SCALING_LOG_10, String.valueOf(cfb.isLog10Scaled())));
        contentList.add(createElement(DimapProductConstants.TAG_NO_DATA_VALUE_USED, String.valueOf(cfb.isNoDataValueUsed())));
        contentList.add(createElement(DimapProductConstants.TAG_NO_DATA_VALUE, String.valueOf(cfb.getNoDataValue())));

        final ArrayList filterBandInfoList = new ArrayList();
        filterBandInfoList.add(createElement(DimapProductConstants.TAG_FILTER_SOURCE, cfb.getSource().getName()));

        final ArrayList filterKernelList = new ArrayList();
        filterKernelList.add(createElement(DimapProductConstants.TAG_KERNEL_WIDTH, String.valueOf(cfb.getKernel().getWidth())));
        filterKernelList.add(createElement(DimapProductConstants.TAG_KERNEL_HEIGHT, String.valueOf(cfb.getKernel().getHeight())));
        filterKernelList.add(createElement(DimapProductConstants.TAG_KERNEL_FACTOR, String.valueOf(cfb.getKernel().getFactor())));
        filterKernelList.add(createElement(DimapProductConstants.TAG_KERNEL_DATA, StringUtils.arrayToCsv(cfb.getKernel().getKernelData(null))));

        final Element filterKernel = new Element(DimapProductConstants.TAG_FILTER_KERNEL);
        filterKernel.addContent(filterKernelList);
        filterBandInfoList.add(filterKernel);

        final Element filterBandInfo = new Element(DimapProductConstants.TAG_FILTER_BAND_INFO);
        filterBandInfo.setAttribute("bandType", "ConvolutionFilterBand");
        filterBandInfo.addContent(filterBandInfoList);
        contentList.add(filterBandInfo);

        final Element root = new Element(DimapProductConstants.TAG_SPECTRAL_BAND_INFO);
        root.setContent(contentList);
        return root;

    }

    private static Element createElement(String tagName, String text) {
        final Element elem = new Element(tagName);
        elem.setText(text);
        return elem;
    }

}
