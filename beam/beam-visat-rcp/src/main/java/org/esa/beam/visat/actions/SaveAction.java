/*
 * $Id: SaveAction.java,v 1.1 2009-04-27 13:08:25 lveci Exp $
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
package org.esa.beam.visat.actions;

import org.esa.beam.framework.dataio.ProductIO;
import org.esa.beam.framework.dataio.ProductReaderPlugIn;
import org.esa.beam.framework.dataio.ProductWriter;
import org.esa.beam.framework.dataio.ProductReader;
import org.esa.beam.framework.datamodel.ProductNode;
import org.esa.beam.framework.ui.command.CommandEvent;
import org.esa.beam.framework.ui.command.ExecCommand;
import org.esa.beam.visat.VisatApp;

/**
 * This action saves the selected product.
 *
 * @author Marco Peters
 * @version $Revision: 1.1 $ $Date: 2009-04-27 13:08:25 $
 */
public class SaveAction extends ExecCommand {

    @Override
    public void actionPerformed(final CommandEvent event) {
        VisatApp.getApp().saveSelectedProduct();
    }

    @Override
    public void updateState(final CommandEvent event) {
        boolean enable = false;
        final ProductNode selectedProductNode = VisatApp.getApp().getSelectedProductNode();
        if (selectedProductNode != null) {
            ProductReader productReader = selectedProductNode.getProductReader();
            if (productReader != null) {
                ProductReaderPlugIn readerPlugIn = productReader.getReaderPlugIn();
                if (readerPlugIn != null) {
                    String[] formatNames = readerPlugIn.getFormatNames();
                    for (int i = 0; i < formatNames.length && !enable; i++) {
                        String formatName = formatNames[i];
                        ProductWriter writer = ProductIO.getProductWriter(formatName);
                        if (writer != null) {
                            enable = true;
                        }
                    }
                } else {
                    // No ReaderPlugIn found so the reader is some kind of AbstractProductBuilder
                    // --> Save should be always anabled
                    enable = true;
                }
            }
        }
        setEnabled(enable);
    }
}
