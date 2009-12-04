package com.bc.ceres.swing.figure;

import com.bc.ceres.grender.Rendering;
import com.bc.ceres.grender.Viewport;

import java.awt.Cursor;
import java.awt.Graphics2D;
import java.awt.Shape;
import java.awt.geom.AffineTransform;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;

public abstract class AbstractHandle extends AbstractFigure implements Handle {

    private final Figure figure;
    private final FigureStyle normalStyle;
    private final FigureStyle selectedStyle;
    private final FigureChangeListener listener;
    private final Point2D.Double location;
    private Shape shape;
    private boolean selected;

    protected AbstractHandle(Figure figure,
                             FigureStyle normalStyle,
                             FigureStyle selectedStyle) {
        this.figure = figure;
        this.normalStyle = normalStyle;
        this.selectedStyle = selectedStyle;

        this.listener = new AbstractFigureChangeListener() {
            @Override
            public void figureChanged(FigureChangeEvent e) {
                updateLocation();
            }
        };
        this.figure.addChangeListener(listener);
        this.location = new Point2D.Double();
    }

    public double getX() {
        return location.x;
    }

    public double getY() {
        return location.y;
    }

    @Override
    public Point2D getLocation() {
        return (Point2D) location.clone();
    }

    public void setLocation(double x, double y) {
        location.setLocation(x, y);
    }

    public abstract void updateLocation();

    public Figure getFigure() {
        return figure;
    }

    public FigureStyle getNormalStyle() {
        return normalStyle;
    }

    @Override
    public Shape getShape() {
        return shape;
    }

    @Override
    public void setShape(Shape shape) {
        this.shape = shape;
    }

    @Override
    public Rank getRank() {
        return Rank.POLYGONAL;
    }

    @Override
    public Rectangle2D getBounds() {
        return shape.getBounds2D();
    }

    /**
     * The default implementation returns {@code true}.
     *
     * @return Always {@code true}.
     */
    @Override
    public boolean isSelectable() {
        return true;
    }

    @Override
    public boolean isSelected() {
        return selected;
    }

    @Override
    public void setSelected(boolean selected) {
        this.selected = selected;
    }

    @Override
    public boolean isCloseTo(Point2D point, AffineTransform m2v) {
        Point2D delta = new Point2D.Double(point.getX() - location.getX(),
                                           point.getY() - location.getY());
        m2v.deltaTransform(delta, delta);
        return getShape().contains(delta);
    }

    @Override
    public void dispose() {
        super.dispose();
        figure.removeChangeListener(listener);
    }

    @Override
    public Cursor getCursor() {
        return Cursor.getPredefinedCursor(Cursor.HAND_CURSOR);
    }

    @Override
    public abstract void move(double dx, double dy);

    @Override
    public final void draw(Rendering rendering) {
        final Graphics2D g = rendering.getGraphics();
        final Viewport vp = rendering.getViewport();
        final AffineTransform oldTransform = g.getTransform();

        try {
            AffineTransform m2v = vp.getModelToViewTransform();
            Point2D transfLocation = m2v.transform(location, null);
            AffineTransform newTransform = new AffineTransform(oldTransform);
            newTransform.concatenate(
                    AffineTransform.getTranslateInstance(transfLocation.getX(), transfLocation.getY()));
            g.setTransform(newTransform);

            drawHandle(g);

        } finally {
            g.setTransform(oldTransform);
        }
    }

    protected void drawHandle(Graphics2D g) {
        FigureStyle handleStyle = isSelected() ? selectedStyle : normalStyle;

        g.setPaint(handleStyle.getFillPaint());
        g.fill(getShape());

        g.setPaint(handleStyle.getStrokePaint());
        g.setStroke(handleStyle.getStroke());
        g.draw(getShape());
    }
}