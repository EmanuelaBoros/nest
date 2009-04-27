/*
 * $Id: ShowMetadataViewAction.java,v 1.1 2009-04-27 13:08:25 lveci Exp $
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
package org.esa.beam.visat.actions;

import org.esa.beam.framework.datamodel.MetadataElement;
import org.esa.beam.framework.datamodel.Product;
import org.esa.beam.framework.datamodel.ProductNode;
import org.esa.beam.framework.datamodel.ProductNodeEvent;
import org.esa.beam.framework.datamodel.ProductNodeListenerAdapter;
import org.esa.beam.framework.ui.UIUtils;
import org.esa.beam.framework.ui.command.CommandEvent;
import org.esa.beam.framework.ui.command.ExecCommand;
import org.esa.beam.framework.ui.product.ProductMetadataView;
import org.esa.beam.framework.ui.product.ProductSceneView;
import org.esa.beam.visat.VisatApp;

import java.awt.Cursor;

import javax.swing.Icon;
import javax.swing.JInternalFrame;

/**
 * This action opens an Metadata View of the currently selected Metadata Node
 *
 * @author Marco Zuehlke
 * @version $Revision: 1.1 $ $Date: 2009-04-27 13:08:25 $
 * @since BEAM 4.6
 */
public class ShowMetadataViewAction extends ExecCommand {
    public static String ID = "showMetadataView";

    @Override
    public void actionPerformed(final CommandEvent event) {
        final VisatApp visatApp = VisatApp.getApp();
        openMetadataView((MetadataElement) visatApp.getSelectedProductNode());
    }
    
    @Override
    public void updateState(final CommandEvent event) {
        setEnabled(VisatApp.getApp().getSelectedProductNode() instanceof MetadataElement);
    }

    public ProductMetadataView openMetadataView(final MetadataElement element) {
        ProductMetadataView metadataView = new ProductMetadataView(element);
        openInternalFrame(metadataView);

        return metadataView;
    }
    
    public JInternalFrame openInternalFrame(ProductMetadataView metadataView) {
        final VisatApp visatApp = VisatApp.getApp();

        visatApp.setStatusBarMessage("Creating metadata view...");
        visatApp.getMainFrame().setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));

        JInternalFrame metadataFrame = null;
        try {
            metadataView.setCommandUIFactory(getCommandUIFactory());
            final Icon icon = UIUtils.loadImageIcon("icons/RsMetaData16.gif");
            final MetadataElement element = metadataView.getMetadataElement();
            metadataFrame = visatApp.createInternalFrame(element.getDisplayName(),
                                                                     icon,
                                                                     metadataView, null);
            final Product product = metadataView.getProduct();
            final JInternalFrame internalFrame = metadataFrame;
            product.addProductNodeListener(new ProductNodeListenerAdapter() {
                @Override
                public void nodeChanged(final ProductNodeEvent event) {
                    if (event.getSourceNode() == element &&
                        event.getPropertyName().equalsIgnoreCase(ProductNode.PROPERTY_NAME_NAME)) {
                        internalFrame.setTitle(element.getDisplayName());
                    }
                }
            });
            updateState();
        } catch (Exception e) {
            visatApp.handleUnknownException(e);
        }

        visatApp.getMainFrame().setCursor(Cursor.getDefaultCursor());
        visatApp.clearStatusBarMessage();

        return metadataFrame;
    }
}
