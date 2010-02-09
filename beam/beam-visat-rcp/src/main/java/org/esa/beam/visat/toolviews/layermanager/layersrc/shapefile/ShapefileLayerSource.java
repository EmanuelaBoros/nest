/*
 * $Id: ShapefileLayerSource.java,v 1.5 2010-02-09 18:03:39 lveci Exp $
 *
 * Copyright (C) 2009 by Brockmann Consult (info@brockmann-consult.de)
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
package org.esa.beam.visat.toolviews.layermanager.layersrc.shapefile;

import org.esa.beam.framework.ui.layer.AbstractLayerSourceAssistantPage;
import org.esa.beam.framework.ui.layer.LayerSource;
import org.esa.beam.framework.ui.layer.LayerSourcePageContext;

/**
 * A layer source for ESRI shape files.
 * <p/>
 * Unstable API. Use at own risk.
 *
 * @author Marco Zuehlke
 * @version $Revision: 1.5 $ $Date: 2010-02-09 18:03:39 $
 * @since BEAM 4.6
 */
public class ShapefileLayerSource implements LayerSource {

    static final String PROPERTY_NAME_FILE_PATH = "fileName";
    static final String PROPERTY_NAME_FEATURE_COLLECTION = "featureCollection";
    static final String PROPERTY_NAME_FEATURE_COLLECTION_CRS = "featureCollectionCrs";
    static final String PROPERTY_NAME_STYLES = "styles";
    static final String PROPERTY_NAME_SELECTED_STYLE = "selectedStyle";

    @Override
    public boolean isApplicable(LayerSourcePageContext pageContext) {
        return pageContext.getAppContext().getSelectedProduct().getGeoCoding() == null;
    }

    @Override
    public boolean hasFirstPage() {
        return true;
    }

    @Override
    public AbstractLayerSourceAssistantPage getFirstPage(LayerSourcePageContext pageContext) {
        return new ShapefileAssistantPage1();
    }

    @Override
    public boolean canFinish(LayerSourcePageContext pageContext) {
        return false;
    }

    @Override
    public boolean performFinish(LayerSourcePageContext pageContext) {
        return false;
    }

    @Override
    public void cancel(LayerSourcePageContext pageContext) {
    }
}
