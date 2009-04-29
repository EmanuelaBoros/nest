/*
 * $Id: DeleteAction.java,v 1.2 2009-04-29 19:35:16 lveci Exp $
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

import org.esa.beam.framework.datamodel.Band;
import org.esa.beam.framework.datamodel.BitmaskDef;
import org.esa.beam.framework.datamodel.Product;
import org.esa.beam.framework.datamodel.ROIDefinition;
import org.esa.beam.framework.datamodel.RasterDataNode;
import org.esa.beam.framework.datamodel.TiePointGrid;
import org.esa.beam.framework.datamodel.VirtualBand;
import org.esa.beam.framework.ui.command.CommandEvent;
import org.esa.beam.framework.ui.command.ExecCommand;
import org.esa.beam.util.Debug;
import org.esa.beam.visat.VisatApp;

import javax.swing.JInternalFrame;
import javax.swing.JOptionPane;
import java.beans.PropertyVetoException;
import java.util.ArrayList;
import java.util.List;

/**
 * This action deletes the currently selected product node, if possible.
 *
 * @author Marco Peters
 * @version $Revision: 1.2 $ $Date: 2009-04-29 19:35:16 $
 */
public class DeleteAction extends ExecCommand {

    @Override
    public void actionPerformed(final CommandEvent event) {
        deleteObject(VisatApp.getApp().getSelectedProductNode());
    }

    @Override
    public void updateState(final CommandEvent event) {
        setEnabled(isDeleteObjectActionPossible());
    }

    private static void deleteObject(final Object object) {
        if (object instanceof Band) {
            final Band band = (Band) object;
            final String[] virtualBands = getVirtualBandsReferencing(band);
            final String[] validMaskNodes = getRasterDataNodesValidMaskReferencing(band);
            final String[] roiNodes = getRasterDataNodesRoiReferencing(band);
            final String[] bitmaskDefs = getBitmaskDefsReferencing(band);
            String message = "Do you really want to delete the band '" + band.getName() + "'?\n"
                             + "This action cannot be undone.\n\n";
            if (virtualBands.length > 0
                || validMaskNodes.length > 0
                || roiNodes.length > 0
                || bitmaskDefs.length > 0) {
                message += "The band to be deleted is referenced by\n"; /*I18N*/
            }
            String indent = "    ";
            if (virtualBands.length > 0) {
                message += "the expression of virtual band(s):\n"; /*I18N*/
                for (String virtualBand : virtualBands) {
                    message += indent + virtualBand + "\n";
                }
            }
            if (validMaskNodes.length > 0) {
                message += "the valid-mask expression of band(s) or tie-point grid(s)\n"; /*I18N*/
                for (String validMaskNode : validMaskNodes) {
                    message += indent + validMaskNode + "\n";
                }
            }
            if (roiNodes.length > 0) {
                message += "the ROI of band(s) or tie-point grid(s):\n"; /*I18N*/
                for (String roiNode : roiNodes) {
                    message += indent + roiNode + "\n";
                }
            }
            if (bitmaskDefs.length > 0) {
                message += "the bitmask definition(s):\n"; /*I18N*/
                for (String bitmaskDef : bitmaskDefs) {
                    message += indent + bitmaskDef + "\n";
                }
            }

            final int status = VisatApp.getApp().showQuestionDialog("Delete Band",
                                                                    message, null);
            if (status == JOptionPane.YES_OPTION) {
                final JInternalFrame[] internalFrames = VisatApp.getApp().findInternalFrames(band);
                for (final JInternalFrame internalFrame : internalFrames) {
                    try {
                        internalFrame.setClosed(true);
                    } catch (PropertyVetoException e) {
                        Debug.trace(e);
                    }
                }
                if (band.hasRasterData()) {
                    band.unloadRasterData();
                }
                final Product product = band.getProduct();
                product.removeBand(band);
            }
        } else {
            VisatApp.getApp().showInfoDialog("Cannot delete the selected object.", null);
        }
    }


    private static boolean isDeleteObjectActionPossible() {
        return VisatApp.getApp().getSelectedProductNode() instanceof Band;
    }

    private static String[] getRasterDataNodesValidMaskReferencing(final RasterDataNode node) {
        final Product product = node.getProduct();
        final List<String> namesList = new ArrayList<String>();
        if (product != null) {
            for (int i = 0; i < product.getNumBands(); i++) {
                final Band band = product.getBandAt(i);
                if (band != node) {
                    if (isNodeReferencedByExpression(node, band.getValidPixelExpression())) {
                        namesList.add(band.getName());
                    }
                }
            }
            for (int i = 0; i < product.getNumTiePointGrids(); i++) {
                final TiePointGrid tiePointGrid = product.getTiePointGridAt(i);
                if (tiePointGrid != node) {
                    if (isNodeReferencedByExpression(node, tiePointGrid.getValidPixelExpression())) {
                        namesList.add(tiePointGrid.getName());
                    }
                }
            }
        }
        return namesList.toArray(new String[namesList.size()]);
    }

    private static String[] getRasterDataNodesRoiReferencing(final RasterDataNode node) {
        final Product product = node.getProduct();
        final List<String> namesList = new ArrayList<String>();
        if (product != null) {
            for (int i = 0; i < product.getNumBands(); i++) {
                final Band band = product.getBandAt(i);
                if (band != node) {
                    final ROIDefinition roiDefinition = band.getROIDefinition();
                    if (roiDefinition != null) {
                        if (isNodeReferencedByExpression(node, roiDefinition.getBitmaskExpr())) {
                            namesList.add(band.getName());
                        }
                    }
                }
            }
            for (int i = 0; i < product.getNumTiePointGrids(); i++) {
                final TiePointGrid tiePointGrid = product.getTiePointGridAt(i);
                if (tiePointGrid != node) {
                    final ROIDefinition roiDefinition = tiePointGrid.getROIDefinition();
                    if (roiDefinition != null) {
                        if (isNodeReferencedByExpression(node, roiDefinition.getBitmaskExpr())) {
                            namesList.add(tiePointGrid.getName());
                        }
                    }
                }
            }
        }
        return namesList.toArray(new String[namesList.size()]);
    }

    private static String[] getBitmaskDefsReferencing(final RasterDataNode node) {
        final Product product = node.getProduct();
        final List<String> namesList = new ArrayList<String>();
        if (product != null) {
            final BitmaskDef[] bitmaskDefs = product.getBitmaskDefs();
            for (final BitmaskDef bitmaskDef : bitmaskDefs) {
                if (isNodeReferencedByExpression(node, bitmaskDef.getExpr())) {
                    namesList.add(bitmaskDef.getName());
                }
            }
        }
        return namesList.toArray(new String[namesList.size()]);
    }

    private static String[] getVirtualBandsReferencing(final RasterDataNode node) {
        final Product product = node.getProduct();
        final List<String> namesList = new ArrayList<String>();
        if (product != null) {
            for (int i = 0; i < product.getNumBands(); i++) {
                final Band band = product.getBandAt(i);
                if (band instanceof VirtualBand) {
                    final VirtualBand virtualBand = (VirtualBand) band;
                    if (isNodeReferencedByExpression(node, virtualBand.getExpression())) {
                        namesList.add(virtualBand.getName());
                    }
                }
            }
        }
        return namesList.toArray(new String[namesList.size()]);
    }

    private static boolean isNodeReferencedByExpression(final RasterDataNode node, final String expression) {
        if (expression == null || expression.trim().length() == 0) {
            return false;
        }
        return expression.matches(".*\\b" + node.getName() + "\\b.*");
    }
}
