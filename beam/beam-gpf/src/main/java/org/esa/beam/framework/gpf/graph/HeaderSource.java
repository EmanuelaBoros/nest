/*
 * $Id: HeaderSource.java,v 1.1 2009-04-28 14:37:14 lveci Exp $
 *
 * Copyright (C) 2008 by Brockmann Consult (info@brockmann-consult.de)
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
package org.esa.beam.framework.gpf.graph;

import com.thoughtworks.xstream.converters.Converter;
import com.thoughtworks.xstream.converters.MarshallingContext;
import com.thoughtworks.xstream.converters.UnmarshallingContext;
import com.thoughtworks.xstream.io.HierarchicalStreamReader;
import com.thoughtworks.xstream.io.HierarchicalStreamWriter;

/**
 * Created by marcoz.
 *
 * @author marcoz
 * @version $Revision: 1.1 $ $Date: 2009-04-28 14:37:14 $
 */
public class HeaderSource {
    private String name;
    private String location;
    private boolean optional = false;
    private String description = "default";
    
    /**
     * Gets the name of source product.
     * @return the name
     */
    public String getName() {
        return name;
    }
    
    /**
     * Sets the name of the source product.
     * 
     * @param name the name to set
     */
    public void setName(String name) {
        this.name = name;
    }

    /**
     * Gets the location of the source product.
     * 
     * @return the source product location
     */
    public String getLocation() {
        return location;
    }

    /**
     * Sets the location of the source product.
     * @param location
     */
    public void setLocation(String location) {
        this.location = location;
    }
    
    /**
     * @return {@code true} if the source product is optional.
     *         In this case the field value thus may be {@code null}.
     *         Defaults to {@code false}.
     */
    public boolean isOptional() {
        return optional;
    }
    
    /**
     * Sets if the source product is optional.

     * @param optional
     */
    public void setOptional(boolean optional) {
        this.optional = optional;
    }
    
    /**
     * Gets the description for this source product.
     * 
     * @return the descriptioon
     */
    public String getDescription() {
        return description;
    }

    /**
     * Sets the description for this source product.
     * 
     * @param description
     */
    public void setDescription(String description) {
        this.description = description;
    }
    
    public static class Converter  implements com.thoughtworks.xstream.converters.Converter {
        
        public boolean canConvert(Class aClass) {
            return HeaderSource.class.equals(aClass);
        }

        @Override
        public void marshal(Object source, HierarchicalStreamWriter writer,
                MarshallingContext context) {
            HeaderSource headerSource = (HeaderSource) source;
            writer.addAttribute("name", headerSource.getName());
            writer.addAttribute("description", headerSource.getDescription());
            writer.addAttribute("optional", Boolean.toString(headerSource.isOptional()));
            writer.setValue(headerSource.getLocation());
            
        }

        @Override
        public Object unmarshal(HierarchicalStreamReader reader,
                UnmarshallingContext context) {
            HeaderSource headerSource = new HeaderSource();
            headerSource.setName(reader.getAttribute("name"));
            headerSource.setDescription(reader.getAttribute("description"));
            headerSource.setOptional(Boolean.parseBoolean(reader.getAttribute("optional")));
            headerSource.setLocation(reader.getValue());
            return headerSource;
        }
    }
}
