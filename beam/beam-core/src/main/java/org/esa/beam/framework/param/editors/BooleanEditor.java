/*
 * $Id: BooleanEditor.java,v 1.1 2009-04-28 14:39:33 lveci Exp $
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
package org.esa.beam.framework.param.editors;

import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;

import javax.swing.JCheckBox;
import javax.swing.JComponent;

import org.esa.beam.framework.param.AbstractParamEditor;
import org.esa.beam.framework.param.Parameter;

/**
 * An editor which uses a {@link javax.swing.JCheckBox}.
 */
public class BooleanEditor extends AbstractParamEditor {

    private JCheckBox _checkBox;

    public BooleanEditor(Parameter parameter) {
        super(parameter, false);
    }

    public JCheckBox getCheckBox() {
        return _checkBox;
    }

    /**
     * Gets the UI component used to edit the parameter's value.
     */
    public JComponent getEditorComponent() {
        return getCheckBox();
    }

    @Override
    protected void initUI() {
        // do not call super.initUI() since we don't want any labels to be created
        _checkBox = new JCheckBox();
        nameEditorComponent(_checkBox);
        if (getParameter().getProperties().getLabel() != null) {
            _checkBox.setText(getParameter().getProperties().getLabel());
        }
        if (getParameter().getProperties().getDescription() != null) {
            _checkBox.setToolTipText(getParameter().getProperties().getDescription());
        }
        _checkBox.addItemListener(new ItemListener() {

            public void itemStateChanged(ItemEvent event) {
                updateParameter();
            }
        });
    }

    @Override
    public void updateUI() {
        super.updateUI();

        boolean newValue;
        if (getParameter().getValue() instanceof Boolean) {
            newValue = (Boolean) getParameter().getValue();
        } else {
            newValue = Boolean.valueOf(getParameter().getValueAsText());
        }
        if (getCheckBox().isSelected() != newValue) {
            getCheckBox().setSelected(newValue);
        }
        if (getCheckBox().isEnabled() != isEnabled()) {
            getCheckBox().setEnabled(isEnabled());
        }
    }

    private void updateParameter() {
        boolean newValue = getCheckBox().isSelected();
        getParameter().setValue(newValue, null);
    }
}
