/*
 * $Id: LayerTypeRegistry.java,v 1.3 2010-02-09 14:10:28 lveci Exp $
 * 
 * Copyright (C) 2009 by Brockmann Consult (info@brockmann-consult.de)
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
package com.bc.ceres.glayer;

import com.bc.ceres.core.ServiceRegistry;
import com.bc.ceres.core.ServiceRegistryListener;
import com.bc.ceres.core.ServiceRegistryManager;
import com.bc.ceres.core.runtime.RuntimeContext;

import java.util.HashMap;
import java.util.Map;
import java.util.ServiceLoader;
import java.util.Set;

/**
 * A registry for layer type instances.
 * <p/>
 * In order to register new layer types, use the standard {@code META-INF/services}
 * JAR service provider interface (SPI). The service provider name is identical to this
 * class' fully qualified name: {@code com.bc.ceres.glayer.LayerType}.
 *
 * @author Marco Zuehlke
 * @since Ceres 0.10
 */
public class LayerTypeRegistry {

    private final ServiceRegistry<LayerType> serviceRegistry;
    private final Map<String, String> aliases;

    /**
     * Return the instance for the given layerType class name.
     *
     * @param name The name of the layer type.
     * @return the instance
     */
    public static LayerType getLayerType(String name) {
        LayerType layerType = getRegistry().getService(name);
        if (layerType != null) {
            return layerType;
        }
        Map<String, String> map = Holder.instance.aliases;
        String layerTypeName = map.get(name);
        if (layerTypeName != null) {
            layerType = getRegistry().getService(layerTypeName);
        }
        return layerType;
    }

    /**
     * Return the instance for the given layerType class.
     *
     * @param layerTypeClass The class of the layer type.
     * @return the instance
     */
    public static <T extends LayerType> T getLayerType(Class<T> layerTypeClass) {
        return (T) getRegistry().getService(layerTypeClass.getName());
    }

    private static ServiceRegistry<LayerType> getRegistry() {
        return Holder.instance.serviceRegistry;
    }

    private LayerTypeRegistry() {
        serviceRegistry = ServiceRegistryManager.getInstance().getServiceRegistry(LayerType.class);
        aliases = new HashMap<String, String>(20);
        serviceRegistry.addListener(new ServiceRegistryListener<LayerType>() {
            public void serviceAdded(ServiceRegistry<LayerType> registry, LayerType layerType) {
                registerAliases(layerType);
            }

            public void serviceRemoved(ServiceRegistry<LayerType> registry, LayerType layerType) {
                unregisterAliases(layerType);
            }
        });
        Set<LayerType> services = serviceRegistry.getServices();
        for (LayerType layerType : services) {
            registerAliases(layerType);
        }
        if (!RuntimeContext.isAvailable()) {
            final ServiceLoader<LayerType> serviceLoader = ServiceLoader.load(LayerType.class);
            for (final LayerType layerType : serviceLoader) {
                serviceRegistry.addService(layerType);
            }
        }
    }
    
    private void registerAliases(LayerType layerType) {
        String layerTypeClassName = layerType.getClass().getName();
        aliases.put(layerType.getName(), layerTypeClassName);
        String[] layerTypeAliases = layerType.getAliases();
        for (String aliasName : layerTypeAliases) {
            aliases.put(aliasName, layerTypeClassName);
        }
    }

    private void unregisterAliases(LayerType layerType) {
        String layerTypeClassName = layerType.getClass().getName();
        String layerTypeName = layerType.getName();
        if (aliases.get(layerTypeName).equalsIgnoreCase(layerTypeClassName)) {
            aliases.remove(layerTypeName);
        }
        String[] keys = aliases.keySet().toArray(new String[0]);
        for (String key : keys) {
            if (aliases.get(key).equalsIgnoreCase(layerTypeClassName)) {
                aliases.remove(key);
            }
        }
    }


    // Initialization on demand holder idiom, see
    // http://en.wikipedia.org/wiki/Initialization_on_demand_holder_idiom
    private static class Holder {
        private static final LayerTypeRegistry instance = new LayerTypeRegistry();
    }
}
