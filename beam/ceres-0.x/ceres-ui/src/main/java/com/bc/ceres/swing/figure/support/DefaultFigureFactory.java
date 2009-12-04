package com.bc.ceres.swing.figure.support;

import com.bc.ceres.swing.figure.FigureFactory;
import com.bc.ceres.swing.figure.FigureStyle;
import com.bc.ceres.swing.figure.PointFigure;
import com.bc.ceres.swing.figure.ShapeFigure;

import java.awt.Shape;
import java.awt.geom.Point2D;

public class DefaultFigureFactory implements FigureFactory {
    @Override
    public PointFigure createPunctualFigure(Point2D geometry, FigureStyle style) {
        return new DefaultPointFigure(geometry, 4.0);
    }

    @Override
    public ShapeFigure createLinealFigure(Shape geometry, FigureStyle style) {
        return new DefaultShapeFigure(geometry, false, style);
    }

    @Override
    public ShapeFigure createPolygonalFigure(Shape geometry, FigureStyle style) {
        return new DefaultShapeFigure(geometry, true, style);
    }
}
