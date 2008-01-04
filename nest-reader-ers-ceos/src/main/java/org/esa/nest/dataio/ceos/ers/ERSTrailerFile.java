package org.esa.nest.dataio.ceos.ers;

import org.esa.nest.dataio.ceos.CeosFileReader;
import org.esa.nest.dataio.ceos.IllegalCeosFormatException;
import org.esa.nest.dataio.ceos.ers.records.ERSTrailerRecord;
import org.esa.nest.dataio.ceos.records.TrailerFileDescriptorRecord;

import javax.imageio.stream.ImageInputStream;
import java.io.IOException;

/*
 * $Id: ERSTrailerFile.java,v 1.1 2008-01-04 16:23:10 lveci Exp $
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

/**
 * * This class represents a trailer file of an Avnir-2 product.
 *
 * @author Marco Peters
 * @version $Revision: 1.1 $ $Date: 2008-01-04 16:23:10 $
 */
class ERSTrailerFile {

    private ERSTrailerRecord _trailerRecord;
    private CeosFileReader _ceosReader;

    public ERSTrailerFile(final ImageInputStream trailerStream) throws IOException,
                                                                          IllegalCeosFormatException {
        _ceosReader = new CeosFileReader(trailerStream);
        // must be created even it is not (yet) used
        // it is needed for positioning the reader correctly
        new TrailerFileDescriptorRecord(_ceosReader);
        _trailerRecord = new ERSTrailerRecord(_ceosReader);
    }

    public int[] getHistogramBinsForBand(final int index) throws IOException,
                                                                 IllegalCeosFormatException {
        return _trailerRecord.getHistogramFor(index);
    }

    public void close() throws IOException {
        _ceosReader.close();
        _ceosReader = null;
        _trailerRecord = null;
    }
}
