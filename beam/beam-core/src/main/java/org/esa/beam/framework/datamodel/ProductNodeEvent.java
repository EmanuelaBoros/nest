/*
 * $Id: ProductNodeEvent.java,v 1.4 2009-12-11 20:46:13 lveci Exp $
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

import java.util.EventObject;

import org.esa.beam.util.Guardian;

/**
 * A product node event informs a product change listener about the source of the notification.
 *
 * @author Norman Fomferra
 * @version $Revision: 1.4 $ $Date: 2009-12-11 20:46:13 $
 * @see ProductNodeListener
 */
public class ProductNodeEvent extends EventObject {

    public final static int NODE_CHANGED = 0;
    public final static int NODE_ADDED = 1;
    public final static int NODE_REMOVED = 2;
    public final static int NODE_DATA_CHANGED = 3;

    private final ProductNodeGroup nodeGroup;
    private final int type;
    private final String propertyName;
    private final Object oldValue;

    /**
     * Constructs a product node event.
     *
     * @param sourceNode The product node on which the Event initially occurred.
     * @param type       the event type
     */
    public ProductNodeEvent(ProductNode sourceNode, int type) {
        this(sourceNode, type, null);
    }

    /**
     * Constructs a product node event.
     *
     * @param sourceNode The product node on which the Event initially occurred.
     * @param type       the event type
     * @param nodeGroup  If event type is NODE_ADDED or NODE_REMOVED this is the parent group to which the source
     *                   node was added to or removed from.
     */
    public ProductNodeEvent(ProductNode sourceNode, int type, ProductNodeGroup nodeGroup) {
        super(sourceNode);
        this.nodeGroup = getNodeGroup(sourceNode, nodeGroup);
        this.propertyName = type == NODE_DATA_CHANGED ? "data" : null;
        this.type = type;
        this.oldValue = null;
    }

    /**
     * Constructs a productEvent object.
     *
     * @param sourceNode   the product node whose property has changed
     * @param propertyName the name of the property that was changed.
     * @param oldValue     the old value.
     */
    public ProductNodeEvent(final ProductNode sourceNode, final String propertyName, final Object oldValue) {
        super(sourceNode);
        Guardian.assertNotNull("propertyName", propertyName);
        this.nodeGroup = getNodeGroup(sourceNode, null);
        this.propertyName = propertyName;
        this.type = NODE_CHANGED;
        this.oldValue = oldValue;
    }

    /**
     * @return The event type.
     */
    public int getType() {
        return type;
    }

    /**
     * @return A reference to the originating product node, i.e. the one who fired the event.
     */
    public final ProductNode getSourceNode() {
        return (ProductNode) getSource();
    }

    /**
     * @return A reference to the group on which a {@link #NODE_ADDED} or {@link #NODE_REMOVED} event occured. May be null.
     */
    public ProductNodeGroup getGroup() {
        return nodeGroup;
    }

    /**
     * Gets the name of the property that was changed.
     *
     * @return the name of the property that was changed.
     */
    public String getPropertyName() {
        return propertyName;
    }

    /**
     * Gets the old value if there is any given in the constructor.
     *
     * @return the old value.
     */
    public Object getOldValue() {
        return oldValue;
    }

    @Override
    public String toString() {
        return String.format("%s [sourceNode=%s, propertyName=%s, type=%d]",
                             getClass().getName(),
                             getSourceNode(),
                             getPropertyName(),
                             getType());
    }


    /**
     * Gets the event type.
     *
     * @return the event type.
     * @deprecated since BEAM 4.7, use {@link #getType()} instead
     */
    @Deprecated
    public final int getId() {
        return type;
    }


    private static ProductNodeGroup getNodeGroup(ProductNode sourceNode, ProductNodeGroup nodeGroup) {
        if (nodeGroup != null) {
            return nodeGroup;
        }
        ProductNode owner = sourceNode.getOwner();
        if (owner instanceof ProductNodeGroup) {
            return (ProductNodeGroup) owner;
        }
        return null;
    }

}
