/*
 * $Id: PointingFactoryRegistry.java,v 1.2 2009-11-04 17:04:32 lveci Exp $
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
package org.esa.beam.framework.datamodel;

import com.bc.ceres.core.ServiceRegistry;
import com.bc.ceres.core.ServiceRegistryManager;

import org.esa.beam.BeamCoreActivator;
import org.esa.beam.util.Guardian;

import java.util.Set;

/**
 * Created by Marco Peters.
 *
 * @author Marco Peters
 * @version $Revision: 1.2 $ $Date: 2009-11-04 17:04:32 $
 */
public class PointingFactoryRegistry {

    private static ServiceRegistry<PointingFactory> typeToFactoryMap;

    private PointingFactoryRegistry() {
        ServiceRegistryManager serviceRegistryManager = ServiceRegistryManager.getInstance();
        typeToFactoryMap = serviceRegistryManager.getServiceRegistry(PointingFactory.class);
        if (!BeamCoreActivator.isStarted()) {
            BeamCoreActivator.loadServices(typeToFactoryMap);
        }
    }

    public static PointingFactoryRegistry getInstance() {
        return Holder.instance;
    }

    public PointingFactory getPointingFactory(String productType) {
        Guardian.assertNotNullOrEmpty("productType", productType);
        Set<PointingFactory> services = typeToFactoryMap.getServices();
        for (PointingFactory descriptor : services) {
            String[] supportedProductTypes = descriptor.getSupportedProductTypes();
            for (String supportedType : supportedProductTypes) {
                if (productType.equalsIgnoreCase(supportedType)) {
                    return descriptor;
                }
            }
        }
        return null;
    }

    public void addFactory(PointingFactory pointingFactory) {
            typeToFactoryMap.addService(pointingFactory);
    }
    
    // Initialization on demand holder idiom
    private static class Holder {
        private static final PointingFactoryRegistry instance = new PointingFactoryRegistry();
    }
}
