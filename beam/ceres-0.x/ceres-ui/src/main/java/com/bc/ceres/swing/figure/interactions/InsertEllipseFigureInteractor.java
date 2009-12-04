package com.bc.ceres.swing.figure.interactions;

import com.bc.ceres.swing.figure.interactions.InsertRectangularFigureInteractor;

import java.awt.geom.Ellipse2D;
import java.awt.geom.RectangularShape;
import java.awt.geom.Point2D;

public class InsertEllipseFigureInteractor extends InsertRectangularFigureInteractor {
    @Override
    protected RectangularShape createRectangularShape(Point2D point) {
        return new Ellipse2D.Double(point.getX(), point.getY(), 0, 0);
    }
}
