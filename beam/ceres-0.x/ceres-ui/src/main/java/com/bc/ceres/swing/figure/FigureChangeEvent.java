/*
 * $Id: FigureChangeEvent.java,v 1.4 2009-12-29 14:18:36 lveci Exp $
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
package com.bc.ceres.swing.figure;

import java.util.EventObject;

/**
 * This event occurs, when a figure has been changed.
 *
 * @author Marco Zuehlke
 * @author Norman Fomferra
 * @since Ceres 0.10
 */
public class FigureChangeEvent extends EventObject {
    /**
     * Possible event types.
     */
    public enum Type {
        /**
         * Figures have been added to a (source-) figure.
         */
        FIGURES_ADDED,
        /**
         * Figures have been removed from a (source-) figure.
         */
        FIGURES_REMOVED,
        /**
         * A (source-) figure has changed.
         */
        FIGURE_CHANGED,
    }

    /**
     * Figures have been added to a (source-) figure.
     */
    public final static Type FIGURES_ADDED = Type.FIGURES_ADDED;
    /**
     * Figures have been removed from a (source-) figure.
     */
    public final static Type FIGURES_REMOVED = Type.FIGURES_REMOVED;
    /**
     * A (source-) figure has changed.
     */
    public final static Type FIGURE_CHANGED = Type.FIGURE_CHANGED;

    private final Type type;
    private final Figure[] figures;

    /**
     * Constructor.
     *
     * @param sourceFigure The source figure which caused the event.
     * @param type         The type of the event.
     * @param figures      The figures added or removed. Should be {@code null} if the event type
     *                     is {@link #FIGURE_CHANGED}.
     */
    public FigureChangeEvent(Figure sourceFigure, Type type, Figure[] figures) {
        super(sourceFigure);
        this.type = type;
        this.figures = figures != null ? figures.clone() : null;
    }

    /**
     * @return The source figure which caused the event.
     */
    public Figure getSourceFigure() {
        return (Figure) getSource();
    }

    /**
     * @return The type of the event.
     */
    public Type getType() {
        return type;
    }

    /**
     * @return The figures added or removed. Returns {@code null} if the event type is {@link #FIGURE_CHANGED}.
     */
    public Figure[] getFigures() {
        return figures;
    }

    @Override
    public String toString() {
        return "FigureChangeEvent [" +
                "source=" + getSourceFigure() +
                ", type=" + type +
                ", #figures=" + (figures == null ? 0 : figures.length) +
                "]";
    }
}
