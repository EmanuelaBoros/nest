/*
 * $Id: MapProjection.java,v 1.2 2010-02-11 17:02:24 lveci Exp $
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
package org.esa.beam.framework.dataop.maptransf;

import org.esa.beam.framework.param.Parameter;
import org.esa.beam.util.Guardian;
import org.esa.beam.util.ObjectUtils;


/**
 * A map projection is a mathematical model for the transformation of locations from a three-dimensional earth surface
 * to a two-dimensional map representation.
 * 
 * @deprecated since BEAM 4.7, use geotools {@link org.geotools.referencing.operation.projection.MapProjection} instead.
 */
@Deprecated
public class MapProjection implements Cloneable {

    private String _name;
    private MapTransform _mapTransform;
    private String _mapUnit;
    private boolean _preDefined;

    public MapProjection(String name, MapTransform mapTransform) {
        this(name, mapTransform, mapTransform.getDescriptor().getMapUnit());
    }

    public MapProjection(String name, MapTransform mapTransform, boolean preDefined) {
        this(name, mapTransform, mapTransform.getDescriptor().getMapUnit(), preDefined);
    }

    public MapProjection(String name, MapTransform mapTransform, String mapUnit) {
        this(name, mapTransform, mapUnit, false);
    }

    public MapProjection(String name, MapTransform mapTransform, String mapUnit, boolean preDefined) {
        Guardian.assertNotNullOrEmpty("name", name);
        Guardian.assertNotNull("mapTransform", mapTransform);
        _name = name;
        _mapTransform = mapTransform;
        _mapUnit = mapUnit;
        _preDefined = preDefined;
    }

    public String getName() {
        return _name;
    }

    public void setName(String name) {
        _name = name;
    }

    public MapTransform getMapTransform() {
        return _mapTransform;
    }

    public void setMapTransform(MapTransform transform) {
        _mapTransform = transform;
    }

    public String getMapUnit() {
        return _mapUnit;
    }

    public void setMapUnit(String mapUnit) {
        _mapUnit = mapUnit;
    }

    public boolean isPreDefined() {
        return _preDefined;
    }

    public void setPreDefined(boolean preDefined) {
        _preDefined = preDefined;
    }

    /**
     * Tests if a user interface is available. The method is a shorthand for
     * <pre>
     *    getMapTransform().getDescriptor().hasTransformUI();
     * </pre>
     *
     * @return <code>true</code> if a user interface is available, in this case the {@link #getMapTransformUI} method
     *         never returns null.
     */
    public boolean hasMapTransformUI() {
        return getMapTransform().getDescriptor().hasTransformUI();
    }

    /**
     * Gets a user interface for editing the transformation properties of this map projection.
     *
     * @return the user interface or null if editing is not supported. The {@link #hasMapTransformUI()} hasTransformUI} method shall return
     *         <code>false</code> in this case.
     */
    public MapTransformUI getMapTransformUI() {
        MapTransform mapTransform = getMapTransform();
        MapTransformDescriptor descriptor = mapTransform.getDescriptor();
        return descriptor.getTransformUI(mapTransform);
    }

    @Override
    public Object clone() {
        try {
            final MapProjection mapProjection = (MapProjection) super.clone();
            mapProjection.setMapTransform(_mapTransform.createDeepClone());
            return mapProjection;
        } catch (CloneNotSupportedException e) {
            throw new IllegalStateException(e);
        }
    }

    /**
     * Alters the underlying map transformation by changing the values of the transform parameters named
     * "semi_major" and "semi_minor" (if any) to the ones of the supplied ellipsoid.
     *
     * @param ellipsoid the ellipsoid
     */
    public void alterMapTransform(Ellipsoid ellipsoid) {
        final MapTransform oldTransform = getMapTransform();
        final Parameter[] parameters = oldTransform.getDescriptor().getParameters();
        final double[] parameterValues = oldTransform.getParameterValues();

        boolean altered = false;
        for (int i = 0; i < parameters.length; i++) {
            if ("semi_minor".equals(parameters[i].getName())) {
                parameterValues[i] = ellipsoid.getSemiMinor();
                altered = true;
            } else if ("semi_major".equals(parameters[i].getName())) {
                parameterValues[i] = ellipsoid.getSemiMajor();
                altered = true;
            }
        }

        if (altered) {
            final MapTransform newTransform = oldTransform.getDescriptor().createTransform(parameterValues);
            setMapTransform(newTransform);
        }
    }

    @Override
    public int hashCode() {
        return getName().hashCode();
    }

    @Override
    public boolean equals(Object obj) {
        if (obj == this) {
            return true;
        }
        if (obj == null) {
            return false;
        }
        if (obj instanceof MapProjection) {
            MapProjection other = (MapProjection) obj;
            return ObjectUtils.equalObjects(other.getName(), getName())
                   && ObjectUtils.equalObjects(other.getMapUnit(), getMapUnit())
                   && ObjectUtils.equalObjects(other.getMapTransform().getDescriptor(),
                                               getMapTransform().getDescriptor())
                   && ObjectUtils.equalObjects(other.getMapTransform().getParameterValues(),
                                               getMapTransform().getParameterValues());
        }
        return false;
    }
}
