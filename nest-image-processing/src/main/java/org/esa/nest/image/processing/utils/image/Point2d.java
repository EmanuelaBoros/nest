package org.esa.nest.image.processing.utils.image;

/**
 * point2d class
 *
 * @author thomas.boudier@snv.jussieu.fr @created 26 aout 2003
 */
public class Point2d {

    /**
     * Description of the Field
     */
    public double x;
    /**
     * Description of the Field
     */
    public double y;

    /**
     * Constructor for the Point2d object
     */
    public Point2d() {
        x = 0.0;
        y = 0.0;
    }

    public Point2d(double x, double y) {
        this.x = x;
        this.y = y;
    }

    /**
     * Description of the Method
     *
     * @param point2d Description of the Parameter
     * @return Description of the Return Value
     */
    public double dist(Point2d point2d) {
        return (Math.sqrt((x - point2d.x) * (x - point2d.x) + (y - point2d.y) * (y - point2d.y)));
    }

    /**
     * Description of the Method
     *
     * @return Description of the Return Value
     */
    @Override
    public String toString() {
        return "x=" + x + " y=" + y;
    }
}
