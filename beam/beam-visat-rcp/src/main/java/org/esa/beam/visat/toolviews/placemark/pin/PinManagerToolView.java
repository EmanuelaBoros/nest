/*
 * $Id: PinManagerToolView.java,v 1.1 2009-04-27 13:08:26 lveci Exp $
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
package org.esa.beam.visat.toolviews.placemark.pin;

import org.esa.beam.framework.datamodel.*;
import org.esa.beam.visat.toolviews.placemark.PlacemarkManagerToolView;
import org.esa.beam.visat.toolviews.placemark.TableModelFactory;

/**
 * A dialog used to manage the list of pins associated with a selected product.
 */
public class PinManagerToolView extends PlacemarkManagerToolView {

    public static final String ID = PinManagerToolView.class.getName();

    public PinManagerToolView() {
        super(PinDescriptor.INSTANCE, new TableModelFactory() {
            public PinTableModel createTableModel(PlacemarkDescriptor placemarkDescriptor, Product product,
                                                  Band[] selectedBands, TiePointGrid[] selectedGrids) {
                return new PinTableModel(placemarkDescriptor, product, selectedBands, selectedGrids);
            }
        });
    }
}
