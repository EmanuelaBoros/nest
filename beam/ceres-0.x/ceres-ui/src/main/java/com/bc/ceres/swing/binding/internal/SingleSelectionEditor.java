/*
 * $Id: SingleSelectionEditor.java,v 1.1 2010-02-10 19:57:11 lveci Exp $
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
package com.bc.ceres.swing.binding.internal;

import com.bc.ceres.binding.PropertyDescriptor;
import com.bc.ceres.binding.ValueSet;
import com.bc.ceres.swing.binding.BindingContext;
import com.bc.ceres.swing.binding.ComponentAdapter;
import com.bc.ceres.swing.binding.PropertyEditor;

import javax.swing.JComboBox;
import javax.swing.JComponent;

/**
 * A one-of-many selection editor using a Combobox.
 *
 * @author Marco Zuehlke
 * @version $Revision: 1.1 $ $Date: 2010-02-10 19:57:11 $
 * @since BEAM 4.6
 */
public class SingleSelectionEditor extends PropertyEditor {

    @Override
    public boolean isValidFor(PropertyDescriptor propertyDescriptor) {
        ValueSet valueSet = propertyDescriptor.getValueSet();
        Class<?> type = propertyDescriptor.getType();
        return valueSet != null && !type.isArray();
    }
    
    @Override
    public JComponent createEditorComponent(PropertyDescriptor propertyDescriptor, BindingContext bindingContext) {
        JComboBox comboBox = new JComboBox();
        ComponentAdapter adapter = new ComboBoxAdapter(comboBox);
        bindingContext.bind(propertyDescriptor.getName(), adapter);
        return comboBox;
    }
}
