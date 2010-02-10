/*
 * $Id: PropertyEditor.java,v 1.1 2010-02-10 19:57:11 lveci Exp $
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
package com.bc.ceres.swing.binding;

import com.bc.ceres.binding.PropertyDescriptor;

import javax.swing.JComponent;
import javax.swing.JLabel;

/**
 * A factory for editors that can edit values of a certain type
 * described by a {@link com.bc.ceres.binding.PropertyDescriptor}.
 *
 * @author Marco Zuehlke
 * @version $Revision: 1.1 $ $Date: 2010-02-10 19:57:11 $
 * @since BEAM 4.6
 */
public abstract class PropertyEditor {
    
    /**
     * Checks if this value editor can be used for the values
     * described by the given property.
     * 
     * @param propertyDescriptor The value descriptor
     * @return true, is this editor can be used for the given value descriptor
     */
    public boolean isValidFor(PropertyDescriptor propertyDescriptor) {
        return false;
    }
    
    /**
     * Creates the editor component for the {@link com.bc.ceres.binding.PropertyDescriptor} and bind it
     * to a {@link com.bc.ceres.binding.PropertyContainer} using the {@link BindingContext}.
     * 
     * @param propertyDescriptor The value descriptor
     * @param bindingContext The binding context
     * @return the editor component
     */
    public JComponent[] createComponents(PropertyDescriptor propertyDescriptor, BindingContext bindingContext) {
        JComponent editorComponent = createEditorComponent(propertyDescriptor, bindingContext);
        
        JLabel label = new JLabel(propertyDescriptor.getDisplayName() + ":");
        Binding binding = bindingContext.getBinding(propertyDescriptor.getName());
        binding.addComponent(label);
        
        return new JComponent[] {editorComponent, label};
    }
    
    /**
     * Creates the editor component together with an optional label for the {@link com.bc.ceres.binding.PropertyDescriptor} and bind it
     * to a {@link com.bc.ceres.binding.PropertyContainer} using the {@link BindingContext}.
     * 
     * If for this editor component a label is applicable, it is return as the second element in the array.
     * 
     * 
     * @param propertyDescriptor The value descriptor
     * @param bindingContext The binding context
     * @return an array containing the editor component as first element and (if applicable) a label as second element
     */
    public abstract JComponent createEditorComponent(PropertyDescriptor propertyDescriptor, BindingContext bindingContext);

}