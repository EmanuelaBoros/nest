/*
 * $Id: InsertPinInteractor.java,v 1.1 2009-12-11 20:46:15 lveci Exp $
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

import org.esa.beam.framework.datamodel.GeoCoding;
import org.esa.beam.framework.datamodel.PinDescriptor;
import org.esa.beam.framework.ui.product.ProductSceneView;
import org.esa.beam.visat.toolviews.placemark.InsertPlacemarkInteractor;

/**
 * A tool used to create (single click), select (single click on a pin) or edit (double click on a pin) the pins
 * displayed in product scene view.
 */
public class InsertPinInteractor extends InsertPlacemarkInteractor {

    public InsertPinInteractor() {
        super(PinDescriptor.INSTANCE);
    }

    @Override
    public void completeInteraction(ProductSceneView sceneView) {
        if (getSelectedPlacemark() != null
                && sceneView.isCurrentPixelPosValid()) {
            GeoCoding geoCoding = sceneView.getRaster().getGeoCoding();
            if (geoCoding != null && geoCoding.canGetGeoPos()) {
                getSelectedPlacemark().setGeoPos(geoCoding.getGeoPos(getSelectedPlacemark().getPixelPos(), null));
            }
        }
    }
}