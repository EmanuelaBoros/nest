/*
 * $Id: ProductWriterPlugIn.java,v 1.1 2009-04-28 14:39:32 lveci Exp $
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
package org.esa.beam.framework.dataio;

/**
 * The <code>ProductWriterPlugIn</code> interface is implemented by data product writer plug-ins.
 * <p/>
 * <p>XMLCoder plug-ins are used to provide meta-information about a particular data format and to create instances of
 * the actual writer objects.
 * <p/>
 * <p> A plug-in can register itself in the <code>ProductIO</code> plug-in registry or it is automatically found during
 * a classpath scan.
 *
 * @author Norman Fomferra
 * @version $Revision: 1.1 $ $Date: 2009-04-28 14:39:32 $
 * @see org.esa.beam.framework.dataio.ProductReaderPlugIn
 */
public interface ProductWriterPlugIn extends ProductIOPlugIn {

    /**
     * Returns an array containing the classes that represent valid output types for this writer.
     * <p/>
     * <p> Intances of the classes returned in this array are valid objects for the <code>setOutput</code> method of the
     * <code>ProductWriter</code> interface (the method will not throw an <code>InvalidArgumentException</code> in this
     * case).
     *
     * @return an array containing valid output types, never <code>null</code>
     *
     * @see ProductWriter#writeProductNodes
     */
    Class[] getOutputTypes();

    /**
     * Creates an instance of the actual product writer class. This method should never return <code>null</code>.
     *
     * @return a new writer instance, never <code>null</code>
     */
    ProductWriter createWriterInstance();

}
