/*
 * Copyright (C) 2010 Array Systems Computing Inc. http://www.array.ca
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation; either version 3 of the License, or (at your option)
 * any later version.
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, see http://www.gnu.org/licenses/
 */
package org.esa.nest.gpf;

import org.esa.beam.framework.datamodel.Band;
import org.esa.beam.framework.datamodel.Product;
import org.esa.beam.framework.datamodel.VirtualBand;
import org.esa.beam.framework.dataop.resamp.ResamplingFactory;
import org.esa.beam.framework.gpf.ui.BaseOperatorUI;
import org.esa.beam.framework.gpf.ui.UIValidation;
import org.esa.beam.framework.ui.AppContext;
import org.esa.nest.datamodel.Unit;
import org.esa.nest.util.DialogUtils;

import javax.swing.*;
import java.awt.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Map;

/**
 * User interface for CreateStackOp
 */
public class CreateStackOpUI extends BaseOperatorUI {

    private final JList mstBandList = new JList();
    private final JList slvBandList = new JList();

    private ArrayList<Integer> defaultMasterBandIndices = new ArrayList<Integer>(2);
    private ArrayList<Integer> defaultSlaveBandIndices = new ArrayList<Integer>(2);

    private final JComboBox resamplingType = new JComboBox(new String[] { "NONE",
                                                                           ResamplingFactory.NEAREST_NEIGHBOUR_NAME,
                                                                           ResamplingFactory.BILINEAR_INTERPOLATION_NAME,
                                                                           ResamplingFactory.CUBIC_CONVOLUTION_NAME });

    private final JComboBox extent = new JComboBox(new String[] {   CreateStackOp.MASTER_EXTENT,
                                                                    CreateStackOp.MIN_EXTENT,
                                                                    CreateStackOp.MAX_EXTENT });

    @Override
    public JComponent CreateOpTab(String operatorName, Map<String, Object> parameterMap, AppContext appContext) {

        initializeOperatorUI(operatorName, parameterMap);
        final JComponent panel = createPanel();
        initParameters();

        return new JScrollPane(panel);
    }

    @Override
    public void initParameters() {

        final String bandNames[] = getBandNames();
        OperatorUIUtils.initBandList(mstBandList, bandNames);
        OperatorUIUtils.initBandList(slvBandList, bandNames);

        OperatorUIUtils.setSelectedListIndices(mstBandList, defaultMasterBandIndices);
        OperatorUIUtils.setSelectedListIndices(slvBandList, defaultSlaveBandIndices);

        resamplingType.setSelectedItem(paramMap.get("resamplingType"));
        extent.setSelectedItem(paramMap.get("extent"));
    }

    @Override
    public UIValidation validateParameters() {

        return new UIValidation(UIValidation.State.OK, "");
    }

    @Override
    public void updateParameters() {

        OperatorUIUtils.updateBandList(mstBandList, paramMap, "masterBandNames");
        OperatorUIUtils.updateBandList(slvBandList, paramMap, "slaveBandNames");

        paramMap.put("resamplingType", resamplingType.getSelectedItem());
        paramMap.put("extent", extent.getSelectedItem());
    }

    private JComponent createPanel() {

        final JPanel contentPane = new JPanel();
        contentPane.setLayout(new GridBagLayout());
        final GridBagConstraints gbc = DialogUtils.createGridBagConstraints();

        contentPane.add(new JLabel("Master Bands:"), gbc);

        gbc.fill = GridBagConstraints.BOTH;
        gbc.gridx = 1;
        contentPane.add(new JScrollPane(mstBandList), gbc);

        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.gridx = 0;
        gbc.gridy++;

        contentPane.add(new JLabel("Slave Bands:"), gbc);
        gbc.fill = GridBagConstraints.BOTH;
        gbc.gridx = 1;
        contentPane.add(new JScrollPane(slvBandList), gbc);
        gbc.gridx = 0;
        gbc.gridy++;

        DialogUtils.addComponent(contentPane, gbc, "Resampling Type:", resamplingType);
        gbc.gridy++;
        DialogUtils.addComponent(contentPane, gbc, "Output Extents:", extent);
        gbc.gridy++;

        DialogUtils.fillPanel(contentPane, gbc);

        return contentPane;
    }

    @Override
    protected String[] getBandNames() {
        final ArrayList<String> bandNames = new ArrayList<String>(5);
        if(sourceProducts != null) {
            boolean masterBandsSelected = false;
            for(Product prod : sourceProducts) {
                if(sourceProducts.length > 1) {

                    final Band[] bands = prod.getBands();
                    for(int i=0; i < bands.length; ++i) {
                        final Band band = bands[i];
                        bandNames.add(band.getName()+"::"+prod.getName());
                        final int index = bandNames.size()-1;

                        if(!(band instanceof VirtualBand)) {

                            if(!masterBandsSelected) {
                                defaultMasterBandIndices.add(index);
                                if(band.getUnit() != null && band.getUnit().equals(Unit.REAL)) {
                                    if(i+1 < bands.length) {
                                        if(bands[i+1].getUnit() != null && bands[i+1].getUnit().equals(Unit.IMAGINARY)) {
                                            defaultMasterBandIndices.add(index+1);
                                        }
                                    }
                                }
                                masterBandsSelected = true;
                            } else if(index > defaultMasterBandIndices.size()) {
                                defaultSlaveBandIndices.add(index);
                            }
                        }
                    }
                } else {
                    bandNames.addAll(Arrays.asList(prod.getBandNames()));
                }
            }
        }
        return bandNames.toArray(new String[bandNames.size()]);
    }

}