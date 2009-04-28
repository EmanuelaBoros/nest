/*
 * $Id: GeoTiffProductReaderPlugIn.java,v 1.1 2009-04-28 14:37:14 lveci Exp $
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
package org.esa.beam.dataio.geotiff;

import com.sun.media.imageioimpl.plugins.tiff.TIFFIFD;
import com.sun.media.imageioimpl.plugins.tiff.TIFFImageMetadata;
import com.sun.media.imageioimpl.plugins.tiff.TIFFImageReader;
import org.esa.beam.framework.dataio.DecodeQualification;
import org.esa.beam.framework.dataio.ProductReader;
import org.esa.beam.framework.dataio.ProductReaderPlugIn;
import org.esa.beam.util.io.BeamFileFilter;

import javax.imageio.ImageIO;
import javax.imageio.ImageReader;
import javax.imageio.stream.ImageInputStream;
import java.io.File;
import java.util.Iterator;
import java.util.Locale;

public class GeoTiffProductReaderPlugIn implements ProductReaderPlugIn {

    private static final String[] FORMAT_NAMES = new String[]{"GeoTIFF"};

    public DecodeQualification getDecodeQualification(Object input) {
        try {
                final File file = Utils.getFile(input);
            final ImageInputStream stream = ImageIO.createImageInputStream(file);

            try {
                return getDecodeQualificationImpl(stream);
            } finally {
                stream.close();
            }
        } catch (Exception ignore) {
            // nothing to do, return value is already UNABLE
        }

        return DecodeQualification.UNABLE;
    }

    static DecodeQualification getDecodeQualificationImpl(ImageInputStream stream) {
        try {
            Iterator<ImageReader> imageReaders = ImageIO.getImageReaders(stream);
            TIFFImageReader imageReader = null;
            while(imageReaders.hasNext()) {
                final ImageReader reader = imageReaders.next();
                if(reader instanceof TIFFImageReader) {
                    imageReader = (TIFFImageReader) reader;
                    break;
                }
            }
            if(imageReader == null)
                return DecodeQualification.UNABLE;
            
            imageReader.setInput(stream);

            final TIFFImageMetadata imageMetadata = (TIFFImageMetadata) imageReader.getImageMetadata(0);
            final TIFFIFD ifd = imageMetadata.getRootIFD();
            final TiffFileInfo info = new TiffFileInfo(ifd);

            if (info.isGeotiff()) {
                return DecodeQualification.INTENDED;
            }
        } catch (Exception ignore) {
            return DecodeQualification.UNABLE;
        }
        return DecodeQualification.SUITABLE;
    }

    public Class[] getInputTypes() {
        return new Class[]{String.class, File.class};
    }

    public ProductReader createReaderInstance() {
        return new GeoTiffProductReader(this);
    }

    public String[] getFormatNames() {
        return FORMAT_NAMES;
    }

    public String[] getDefaultFileExtensions() {
        return new String[]{".tif", ".tiff"};
    }

    public String getDescription(Locale locale) {
        return "GeoTIFF data product.";
    }

    public BeamFileFilter getProductFileFilter() {
        return new BeamFileFilter(FORMAT_NAMES[0], getDefaultFileExtensions(), getDescription(null));
    }
}
