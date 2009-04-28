/*
 * $Id: OperatorConfiguration.java,v 1.1 2009-04-28 14:37:14 lveci Exp $
 * 
 * Copyright (C) 2008 by Brockmann Consult (info@brockmann-consult.de)
 * 
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation. This program is distributed in the hope it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 * 
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place - Suite 330, Boston, MA 02111-1307, USA.
 */
package org.esa.beam.framework.gpf.internal;

import com.bc.ceres.binding.dom.DomElement;
import org.esa.beam.framework.gpf.Operator;

import java.util.Set;

/**
 * Created by marcoz.
 *
 * @author marcoz
 * @version $Revision: 1.1 $ $Date: 2009-04-28 14:37:14 $
 */
public class OperatorConfiguration {

    private final DomElement configuration;
    private Set<Reference> referenceSet;

    public OperatorConfiguration(DomElement configuration,
                                 Set<Reference> references) {
        this.configuration = configuration;
        this.referenceSet = references;
    }

    public DomElement getConfiguration() {
        return configuration;
    }

    public Set<Reference> getReferenceSet() {
        return referenceSet;
    }

    public static interface Reference {
        public Object getValue();

        public String getParameterName();
    }

    public static class PropertyReference implements Reference {
        final String parameterName;
        final String propertyName;
        final Operator operator;

        public PropertyReference(String parameterName, String propertyName,
                                 Operator operator) {
            this.parameterName = parameterName;
            this.propertyName = propertyName;
            this.operator = operator;
        }

        @Override
        public Object getValue() {
            return operator.getTargetProperty(propertyName);
        }

        @Override
        public String getParameterName() {
            return parameterName;
        }
    }

    public static class ParameterReference implements Reference {

        private final String name;
        private final Object value;

        public ParameterReference(String name, Object value) {
            this.name = name;
            this.value = value;
        }

        @Override
        public String getParameterName() {
            return name;
        }

        @Override
        public Object getValue() {
            return value;
        }

    }
}
