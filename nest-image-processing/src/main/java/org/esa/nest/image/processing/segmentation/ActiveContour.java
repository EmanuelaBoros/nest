/*
 * Copyright (C) 2012 by Array Systems Computing Inc. http://www.array.ca
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation; either version 3 of the License, or (at your option)
 * any later version.
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, see http://www.gnu.org/licenses/
 */
package org.esa.nest.image.processing.segmentation;

import ij.*;
import ij.gui.*;
import ij.process.*;
import java.awt.*;
import java.io.*;
import java.text.NumberFormat;
import java.util.Locale;
import org.esa.nest.image.processing.utils.image.Point2d;
import org.esa.nest.image.processing.utils.configuration.ActiveContourConfiguration;

/**
 * Active Contour finds a contour that best approximates the perimeter of an
 * object
 *
 * @author thomas.boudier@snv.jussieu.fr @created August 26th 2003
 */
public class ActiveContour {

    private Point2d points[];
    private Point2d normale[];
    private Point2d deplace[];
    private double dataDistance;
    private double lambda[];
    private int state[];
    private int NPT;
    private int NMAX = 50000;
    private int OFF;
    private int block;
    private int elimination;
    private boolean closed;
    private ActiveContourConfiguration configuration;
    private ImageProcessor gradientImage;
    private ImageProcessor originalImage;

    /**
     * Constructor
     */
    public ActiveContour() {
    }

    /**
     * Description of the Method
     */
    public void kill() {
        points = null;
        normale = null;
        deplace = null;
        lambda = null;
        state = null;
        System.gc();
    }

    /**
     * Sets the config attribute of the SnakeD object
     *
     * @param sc The new config value
     */
    public void setConfig(ActiveContourConfiguration sc) {
        configuration = sc;
    }

    /**
     * Get number of points
     *
     * @return The nbPoints value
     */
    public int getNbPoints() {
        return NPT;
    }

    /**
     * Gets the point attribute of the ActiveContour object
     *
     * @param i Description of the Parameter
     * @return The point value
     */
    public Point2d getPoint(int i) {
        return points[i];
    }

    /**
     * Gets the points attribute of the ActiveContour object
     *
     * @return The points value
     */
    public Point2d[] getPoints() {
        return points;
    }

    /**
     * Gets the configuration attribute of the ActiveContour object
     *
     * @return The configuration value
     */
    public ActiveContourConfiguration getConfiguration() {
        return configuration;
    }

    /**
     * Gets the lambda attribute of the ActiveContour object
     *
     * @return The lambda value
     */
    public double[] getLambda() {
        return lambda;
    }

    /**
     * Gets the displacement attribute of the ActiveContour object
     *
     * @return The displacement value
     */
    public Point2d[] getDisplacement() {
        return deplace;
    }

    /**
     * Is the snake isClosed
     *
     * @return Description of the Return Value
     */
    public boolean isClosed() {
        return closed;
    }

    public void setOriginalImage(ImageProcessor originalImage) {
        this.originalImage = originalImage;
    }

    /**
     * Draw the snake
     *
     * @param A Description of the Parameter
     * @param col Description of the Parameter
     * @param linewidth Description of the Parameter
     */
    public void drawContour(ImageProcessor A, Color col, int linewidth) {
        int i;
        int x;
        int y;
        int xx;
        int yy;
        A.setColor(col);
        A.setLineWidth(linewidth);
        for (i = 0; i < NPT - 1; i++) {
            x = (int) (points[i].x);
            y = (int) (points[i].y);
            xx = (int) (points[i + 1].x);
            yy = (int) (points[i + 1].y);
            A.drawLine(x, y, xx, yy);
        }
        if (this.isClosed()) {
            x = (int) (points[NPT - 1].x);
            y = (int) (points[NPT - 1].y);
            xx = (int) (points[0].x);
            yy = (int) (points[0].y);
            A.drawLine(x, y, xx, yy);
        }
    }

    /**
     * write output in the FreeD format
     *
     * @param nb section number
     */
    public void writeInFreeD(int nb) {
        try {
            File fichier = new File("freed" + (nb + 1) + ".txt");
            FileWriter fw = new FileWriter(fichier);
            BufferedWriter bw = new BufferedWriter(fw);
            bw.write("MODELITEM cervo \"section " + (nb + 1) + "\";");
            bw.write("\nMODELITEMDATA  \n");
            for (int i = 0; i < NPT; i++) {
                bw.write("" + (int) points[i].x + "," + (int) points[i].y + "  \n");
            }
            bw.write(";\n");
            bw.close();
            fw.close();
        } catch (IOException e) {
        }
    }

    /**
     * write output in a text format
     *
     * @param nom file name
     * @param nb slice number
     */
    public void writeCoordinates(String nom, int nb, double resXY) {
        NumberFormat nf = NumberFormat.getInstance(Locale.ENGLISH);
        nf.setMaximumFractionDigits(3);
        try {
            File fichier = new File(nom + nb + ".txt");
            FileWriter fw = new FileWriter(fichier);
            BufferedWriter bw = new BufferedWriter(fw);
            bw.write("nb\tX\tY\tZ\tXcal\tYcal\n");
            for (int i = 0; i < NPT; i++) {
                bw.write(i + "\t" + nf.format(points[i].x) + "\t" + nf.format(points[i].y) + "\t" + nb + "\t" + nf.format(points[i].x * resXY) + "\t" + nf.format(points[i].y) + "\n");
            }
            bw.close();
            fw.close();
        } catch (IOException e) {
        }
    }

    /**
     * creation of a polygon ROI
     *
     * @param imp image plus
     * @return roi
     */
    PolygonRoi createROI() {
        int xx[] = new int[NPT];
        int yy[] = new int[NPT];
        for (int i = 0; i < NPT; i++) {
            xx[i] = (int) (points[i].x);
            yy[i] = (int) (points[i].y);
        }
        return new PolygonRoi(xx, yy, NPT - 1, Roi.FREEROI);
    }

    /**
     * Initialization of the snake points
     *
     * @param R ROI
     */
    public void initActiveContour(Roi roi) {
        Double pos;
        double Rx;
        double Ry;
        int i = 1;
        double a;
        NPT = 0;

        points = new Point2d[NMAX];
        normale = new Point2d[NMAX];
        deplace = new Point2d[NMAX];
        dataDistance = 0.0;
        state = new int[NMAX];
        lambda = new double[NMAX];

        for (i = 0; i < NMAX; i++) {
            points[i] = new Point2d();
            normale[i] = new Point2d();
            deplace[i] = new Point2d();
        }


        //Calcul des points de la ROI
        if ((roi.getType() == Roi.OVAL) || (roi.getType() == Roi.RECTANGLE)) {
            closed = true;
            Rectangle Rect = roi.getBounds();
            int xc = Rect.x + Rect.width / 2;
            int yc = Rect.y + Rect.height / 2;
            Rx = ((double) Rect.width) / 2;
            Ry = ((double) Rect.height) / 2;
            double theta = 4.0 / (Rx + Ry);
            i = 0;
            for (a = 2 * Math.PI; a > 0; a -= theta) {
                points[i].x = (int) (xc + Rx * Math.cos(a));
                points[i].y = (int) (yc + Ry * Math.sin(a));
                state[i] = 0;
                i++;
            }
            NPT = i;
        } else if (roi.getType() == Roi.LINE) {
            closed = false;
            Line l = (Line) (roi);
            Rx = (l.x2 - l.x1);
            Ry = (l.y2 - l.y1);
            a = Math.sqrt(Rx * Rx + Ry * Ry);
            Rx /= a;
            Ry /= a;
            int ind = 0;
            for (i = 0; i <= l.getLength(); i++) {
                points[ind].x = l.x1 + Rx * i;
                points[ind].y = l.y1 + Ry * i;
                state[ind] = 0;
                ind++;
            }
            NPT = ind;
        } else if ((roi.getType() == Roi.FREEROI) || (roi.getType() == Roi.POLYGON)) {
            closed = true;
            PolygonRoi p = (PolygonRoi) (roi);
            Rectangle rectBound = p.getBounds();
            int NBPT = p.getNCoordinates();
            int pointsX[] = p.getXCoordinates();
            int pointsY[] = p.getYCoordinates();
            for (i = 0; i < NBPT; i++) {
                points[i].x = pointsX[i] + rectBound.x;
                points[i].y = pointsY[i] + rectBound.y;

            }
            NPT = NBPT;
            if (roi.getType() == Roi.POLYGON) {
                this.resample(true);
            }
        } else {
            IJ.showStatus("Selection type not supported");
        }
        block = 0;
        elimination = 0;
        OFF = 0;
    }

    /**
     * regularisation of distance between points
     */
    void resample(boolean init) {
        Point2d temp[];
        Point2d Ta;
        int i;
        int j;
        int k;
        int ii;
        int aj;
        double Dmoy;
        double Dtot;
        double DD;
        double Dmin;
        double Dmax;
        double Di;
        double Dmoyg;
        double normtan;
        double D;
        double D1;

        temp = new Point2d[NMAX];
        Ta = new Point2d();

        Dtot = 0.0;
        Dmin = 1000.0;
        Dmax = 0.0;
        for (i = 1; i < NPT; i++) {
            Di = distance(i, i - 1);
            Dtot += Di;
            if (Di < Dmin) {
                Dmin = Di;
            }
            if (Di > Dmax) {
                Dmax = Di;
            }
        }
        if (((Dmax / Dmin) > 3.0) || (init)) {
            Dmoyg = 1.0;
            temp[0] = new Point2d();
            temp[0].x = points[0].x;
            temp[0].y = points[0].y;
            i = 1;
            ii = 1;
            temp[ii] = new Point2d();
            while (i < NPT) {
                Dmoy = Dmoyg;
                DD = distance(i, i - 1);
                if (DD > Dmoy) {
                    aj = (int) (DD / Dmoy);
                    Ta.x = points[i].x - points[i - 1].x;
                    Ta.y = points[i].y - points[i - 1].y;
                    normtan = Math.sqrt(Ta.x * Ta.x + Ta.y * Ta.y);
                    Ta.x /= normtan;
                    Ta.y /= normtan;
                    for (k = 1; k <= aj; k++) {
                        temp[ii].x = points[i - 1].x + k * Dmoy * Ta.x;
                        temp[ii].y = points[i - 1].y + k * Dmoy * Ta.y;
                        ii++;
                        temp[ii] = new Point2d();
                    }
                }
                i++;
                if ((DD <= Dmoy) && (i < NPT - 1)) {
                    j = i - 1;
                    D = 0.0;
                    while ((D < Dmoy) && (j < NPT - 1)) {
                        D += distance(j, j + 1);
                        j++;
                    }
                    temp[ii].x = points[j].x;
                    temp[ii].y = points[j].y;
                    ii++;
                    temp[ii] = new Point2d();
                    i = j + 1;
                }
                if (i == NPT - 1) {
                    i = NPT;
                }
            }
            temp[ii].x = points[NPT - 1].x;
            temp[ii].y = points[NPT - 1].y;
            NPT = ii + 1;
            for (i = 0; i < NPT; i++) {
                points[i].x = temp[i].x;
                points[i].y = temp[i].y;
            }
        }
    }

    /**
     * main calculus function (matrix inversion)
     *
     * @param iFirstRow first row
     * @param iLastRow last row
     */
    public void calculus(int iFirstRow, int iLastRow) {
        int i;
        Point2d bi;
        Point2d temp;
        Point2d debtemp;
        double mi;
        double gi;
        double di;
        double omega;

        omega = 1.8;
        bi = new Point2d();
        temp = new Point2d();
        debtemp = new Point2d();

        debtemp.x = points[iFirstRow].x;
        debtemp.y = points[iFirstRow].y;

        for (i = iFirstRow; i < iLastRow; i++) {
            bi.x = points[i].x + deplace[i].x;
            bi.y = points[i].y + deplace[i].y;
            //gi = -lambda[i] * lambda[i + 1] - (lambda[i] * lambda[i]);
            //di = -lambda[i] * lambda[i + 1] - (lambda[i + 1] * lambda[i + 1]);
            //mi = (lambda[i] * lambda[i]) + 2.0 * lambda[i] * lambda[i + 1] + (lambda[i + 1] * lambda[i + 1]) + 1.0;
            gi = -lambda[i];
            di = -lambda[i + 1];
            mi = lambda[i] + lambda[i + 1] + 1.0;
            if (i > iFirstRow) {
                temp.x = mi * points[i].x + omega * (-gi * points[i - 1].x - mi * points[i].x - di * points[i + 1].x + bi.x);
                temp.y = mi * points[i].y + omega * (-gi * points[i - 1].y - mi * points[i].y - di * points[i + 1].y + bi.y);
            }
            if ((i == iFirstRow) && (closed)) {
                temp.x = mi * points[i].x + omega * (-gi * points[iLastRow].x - mi * points[i].x - di * points[i + 1].x + bi.x);
                temp.y = mi * points[i].y + omega * (-gi * points[iLastRow].y - mi * points[i].y - di * points[i + 1].y + bi.y);
            }
            if ((i == iFirstRow) && (!closed)) {
                temp.x = points[iFirstRow].x * mi;
                temp.y = points[iFirstRow].y * mi;
            }
            points[i].x = temp.x / mi;
            points[i].y = temp.y / mi;
        }
        // LAST POINT
        if (closed) {
            i = iLastRow;
            bi.x = points[i].x + deplace[i].x;
            bi.y = points[i].y + deplace[i].y;
            //gi = -lambda[i] * lambda[deb] - (lambda[i] * lambda[i]);
            //di = -lambda[i] * lambda[deb] - (lambda[deb] * lambda[deb]);
            //mi = (lambda[i] * lambda[i]) + 2.0 * lambda[i] * lambda[deb] + (lambda[deb] * lambda[deb]) + 1.0;
            gi = -lambda[i];
            di = -lambda[iFirstRow];
            mi = lambda[i] + lambda[iFirstRow] + 1.0;
            temp.x = mi * points[i].x + omega * (-gi * points[i - 1].x - mi * points[i].x - di * debtemp.x + bi.x);
            temp.y = mi * points[i].y + omega * (-gi * points[i - 1].y - mi * points[i].y - di * debtemp.y + bi.y);
            points[i].x = temp.x / mi;
            points[i].y = temp.y / mi;
        }
    }

    /**
     * Description of the Method
     *
     * @return Description of the Return Value
     */
    public double computeDisplacements() {
        
        double sum = 0.0;
        double threshold =  configuration.getGradThreshold();
        double DivForce = configuration.getMaxDisplacement();
        Point2d displ = new Point2d();
        double force;
        sum = 0;
        for (int i = 0; i < NPT; i++) {
            displ.x = 0.0;
            displ.y = 0.0;
            displ = searchTheClosestEdge(i, threshold , 1000, 1000, 0);

            force = Math.sqrt(displ.x * displ.x + displ.y * displ.y);
            if (force > DivForce) {
                deplace[i].x = DivForce * (displ.x / force);
                deplace[i].y = DivForce * (displ.y / force);
            } else {
                deplace[i].x = displ.x;
                deplace[i].y = displ.y;
            }
            force = Math.sqrt(deplace[i].x * deplace[i].x + deplace[i].y * deplace[i].y);

            sum += force;
        }
        return sum;
    }

    /**
     * Description of the Method
     *
     * @param image Description of the Parameter
     */
    public void computeGrad(ImageProcessor image) {
        gradientImage = grad2d_deriche(image, configuration.getAlpha());
    }

    /**
     * search for the closest edge along the normal direction
     *
     * @param iContourPoints number for the snake point
     * @param dEdgeThreshold threshold
     * @param directions directions to look for
     * @return the displacement vector towards the edges
     */
    Point2d searchTheClosestEdge(int iContourPoints, double dEdgeThreshold,
            double dDistancePlus, double dDistanceMinus, int dir) {
        double iy;
        double ix;
        double deplus;
        double demoins;
        double scaleint = configuration.getMaxSearch();
        double Dist;
        double crp = Double.NaN;
        double crm = Double.NaN;
        double bres;
        double ares;
        double bden;
        double bnum;
        double ii;
        Point2d displacement;
        Point2d pos;
        Point2d norm;
        int scale = 10;
        double image_line[] = new double[(int) (2 * scale * scaleint + 1)];

        pos = points[iContourPoints];
        norm = normale[iContourPoints];

        displacement = new Point2d();
        //recherche des points de la normale au point de contour
        int index = 0;
        double step = 1.0 / (double) scale;
        double deb = -scaleint;
        for (ii = deb; ii < scaleint; ii += step) {
            iy = (pos.y + norm.y * ii);
            ix = (pos.x + norm.x * ii);
            if (ix < 0) {
                ix = 0;
            }
            if (iy < 0) {
                iy = 0;
            }
            if (ix >= gradientImage.getWidth()) {
                ix = gradientImage.getWidth() - 1;
            }
            if (iy >= gradientImage.getHeight()) {
                iy = gradientImage.getHeight() - 1;
            }
            image_line[index] = gradientImage.getInterpolatedPixel(ix, iy);
            index++;
        }

        // polygon crossing, avoid self-intersecting snake
        for (int i = 0; i < NPT - 1; i++) {
            if ((i != iContourPoints) && (i != iContourPoints - 1)) {
                bden = (-norm.x * points[i + 1].y + norm.x * points[i].y + norm.y * points[i + 1].x - norm.y * points[i].x);
                bnum = (-norm.x * pos.y + norm.x * points[i].y + norm.y * pos.x - norm.y * points[i].x);
                if (bden != 0) {
                    bres = (bnum / bden);
                } else {
                    bres = 5.0;
                }
                if ((bres >= 0.0) && (bres <= 1.0)) {
                    ares = (float) (-(-points[i + 1].y * pos.x + points[i + 1].y * points[i].x + points[i].y * pos.x + pos.y * points[i + 1].x - pos.y * points[i].x - points[i].y * points[i + 1].x) / (-norm.x * points[i + 1].y + norm.x * points[i].y + norm.y * points[i + 1].x - norm.y * points[i].x));
                    if ((ares > 0.0) && (ares < crp)) {
                        crp = ares;
                    }
                    if ((ares < 0.0) && (ares > crm)) {
                        crm = ares;
                    }
                }
            }
        }
        double coeff_crossing = 0.9;
        crp = crp * coeff_crossing;
        crm = crm * coeff_crossing;

        deplus = Double.POSITIVE_INFINITY;
        demoins = Double.NEGATIVE_INFINITY;

        boolean edge_found = false;
        for (index = 1; index < 2 * scale * scaleint - 1; index++) {
            // check edge threshold
            // local maximum
            if ((image_line[index] >= dEdgeThreshold) && (image_line[index] >= image_line[index - 1]) && (image_line[index] >= image_line[index + 1])) {
                Dist = index * step + deb;
                if ((Dist < 0) && (Dist > demoins)) {
                    demoins = Dist;
                    edge_found = true;
                }
                if ((Dist >= 0) && (Dist < deplus)) {
                    deplus = Dist;
                    edge_found = true;
                }
            }
        }
        state[iContourPoints] = 0;
        //posplus = deplus;
        //posmoins = demoins;

        // no edge found
        if (!edge_found) {
            displacement.x = 0.0;
            displacement.y = 0.0;

            return displacement;
        }

        // check edges found against threshold distances plus and minus
        if (deplus > dDistancePlus) {
            deplus = 2 * scaleint;
        }
        if (demoins < -dDistanceMinus) {
            demoins = -2 * scaleint;
        }
        if (Double.isInfinite(deplus) && Double.isInfinite(demoins)) {
            displacement.x = 0.0;
            displacement.y = 0.0;

            return displacement;
        }

        //System.out.println("num =" + num + " " + demoins + " " + dist_minus + " " + deplus + " " + dist_plus);

        int direction;
        // go to closest edge
        if (-demoins < deplus) {
            displacement.x = norm.x * demoins;
            displacement.y = norm.y * demoins;
            direction = -1;
        } else {
            displacement.x = norm.x * deplus;
            displacement.y = norm.y * deplus;
            direction = 1;
        }

        // test direction
        /*
         * if (((dir != 0) && (dir !=
         * direction))||(originalImage.getPixel((int)pos.x,(int)pos.y)==0)) {
         * displacement.x = 0.0; displacement.y = 0.0;
         *
         * return displacement; }
         *
         */



        /*
         * if (Math.abs(deplus) < Math.abs(demoins)) { displacement.x = norm.x *
         * (double) (deplus); displacement.y = norm.y * (double) (deplus); }
         * else { displacement.x = norm.x * (double) (demoins); displacement.y =
         * norm.y * (double) (demoins); }
         *
         *
         *
         * //snake displacement // if no edge, or crossing, then no
         * displacement if (((deplus == 1000) && (demoins == -1000)) ||
         * ((posplus > crp) && (posmoins < crm)) || ((deplus == 1000) &&
         * (posmoins < crm)) || ((demoins == -1000) && (posplus > crp))) {
         * displacement.x = 0.0; displacement.y = 0.0; } // positive
         * displacement if (((deplus < 1000) && (posplus < crp)) && ((demoins ==
         * -1000) || (posmoins < crm) || (posplus < -posmoins))) {
         * displacement.x = norm.x * (double) (deplus); displacement.y = norm.y
         * * (double) (deplus); } // negative displacement if (((demoins >
         * -1000) && (posmoins > crm)) && ((deplus == 1000) || (posplus > crp)
         * || (posmoins > -posplus))) { displacement.x = norm.x * (double)
         * (demoins); displacement.y = norm.y * (double) (demoins); }
         */


        return displacement;
    }

    /**
     * Description of the Method
     */
    public void compute_normales() {
        for (int i = 0; i < NPT; i++) {
            normale[i] = compute_normale(i);
        }
    }

    /**
     * Description of the Method
     */
    public void compute_lambdas() {
        double force;
        double maxforce = 0.0;
        double minr = configuration.getRegMin();
        double maxr = configuration.getRegMax();

        for (int i = 0; i < NPT; i++) {
            force = Math.sqrt(deplace[i].x * deplace[i].x + deplace[i].y * deplace[i].y);
            if (force > maxforce) {
                maxforce = force;
            }
        }

        for (int i = 0; i < NPT; i++) {
            force = Math.sqrt(deplace[i].x * deplace[i].x + deplace[i].y * deplace[i].y);
            lambda[i] = maxr / (1.0 + ((maxr - minr) / minr) * (force / maxforce));
        }
    }

    /**
     * compute normale
     *
     * @param np number for the snake point
     * @return normal vector
     */
    Point2d compute_normale(int np) {
        Point2d norma;
        Point2d tan;
        double normtan;

        tan = new Point2d();
        norma = new Point2d();

        if (np == 0) {
            if (closed) {
                tan.x = points[1].x - points[NPT - 1].x;
                tan.y = points[1].y - points[NPT - 1].y;
            } else {
                tan.x = points[1].x - points[0].x;
                tan.y = points[1].y - points[0].y;
            }
        }
        if (np == NPT - 1) {
            if (closed) {
                tan.x = points[0].x - points[NPT - 2].x;
                tan.y = points[0].y - points[NPT - 2].y;
            } else {
                tan.x = points[NPT - 1].x - points[NPT - 2].x;
                tan.y = points[NPT - 1].y - points[NPT - 2].y;
            }
        }
        if ((np > 0) && (np < NPT - 1)) {
            tan.x = points[np + 1].x - points[np - 1].x;
            tan.y = points[np + 1].y - points[np - 1].y;
        }
        normtan = Math.sqrt(tan.x * tan.x + tan.y * tan.y);
        if (normtan > 0.0) {
            tan.x /= normtan;
            tan.y /= normtan;
            norma.x = -tan.y;
            norma.y = tan.x;
        }
        return (norma);
    }

    /**
     * destruction
     */
    void destroysnake() {
        Point2d temp[];
        Point2d fo[];
        double lan[];
        int state[];
        int i;
        int j;

        temp = new Point2d[NPT];
        fo = new Point2d[NPT];
        lan = new double[NPT];
        state = new int[NPT];

        j = 0;
        for (i = 0; i < NPT; i++) {
            if (state[i] != 1) {
                temp[j] = new Point2d();
                temp[j].x = points[i].x;
                temp[j].y = points[i].y;
                state[j] = state[i];
                fo[j] = new Point2d();
                fo[j].x = deplace[i].x;
                fo[j].y = deplace[i].y;
                lan[j] = lambda[i];
                j++;
            }
        }
        NPT = j;

        for (i = 0; i < NPT; i++) {
            points[i].x = temp[i].x;
            points[i].y = temp[i].y;
            state[i] = state[i];
            deplace[i].x = fo[i].x;
            deplace[i].y = fo[i].y;
            lambda[i] = lan[i];
        }
    }

    /**
     * distance between two points of the snake
     *
     * @param a number of first point
     * @param b number of second point
     * @return distance
     */
    double distance(int a, int b) {
        return (Math.sqrt(Math.pow(points[a].x - points[b].x, 2.0) + Math.pow(points[a].y - points[b].y, 2.0)));
    }

    /**
     * compute new positions of the snake
     */
    void new_positions() {
        calculus(0, NPT - 1);
    }

    /**
     * Deriche filtering
     *
     * @param iDep image
     * @param alphaD Description of the Parameter
     * @return Description of the Return Value
     */
    private ImageProcessor grad2d_deriche(ImageProcessor iDep, double alphaD) {
        //ImageProcessor iGrad = iDep.createProcessor(iDep.getWidth(), iDep.getHeight());
        ImageProcessor iGrad = new ByteProcessor(iDep.getWidth(), iDep.getHeight());
        int nmem;
        float[] nf_grx = null;
        float[] nf_gry = null;
        float[] a1 = null;
        float[] a2 = null;
        float[] a3 = null;
        float[] a4 = null;
        byte[] result_array = null;

        int icolonnes = 0;
        int icoll = 0;
        int lignes;
        int colonnes;
        int lig_1;
        int lig_2;
        int lig_3;
        int col_1;
        int col_2;
        int col_3;
        int jp1;
        int jm1;
        int im1;
        int ip1;
        int icol_1;
        int icol_2;
        int i;
        int j;
        float ad1;
        float ad2;
        float wd;
        float gzr;
        float gun;
        float an1;
        float an2;
        float an3;
        float an4;
        float an11;

        lignes = iDep.getHeight();
        colonnes = iDep.getWidth();
        nmem = lignes * colonnes;

        lig_1 = lignes - 1;
        lig_2 = lignes - 2;
        lig_3 = lignes - 3;
        col_1 = colonnes - 1;
        col_2 = colonnes - 2;
        col_3 = colonnes - 3;

        /*
         * alloc temporary buffers
         */
        nf_grx = new float[nmem];
        nf_gry = new float[nmem];

        a1 = new float[nmem];
        a2 = new float[nmem];
        a3 = new float[nmem];
        a4 = new float[nmem];

        ad1 = (float) -Math.exp(-alphaD);
        ad2 = 0;
        an1 = 1;
        an2 = 0;
        an3 = (float) Math.exp(-alphaD);
        an4 = 0;
        an11 = 1;

        /*
         * FIRST STEP: Y GRADIENT
         */
        /*
         * x-smoothing
         */
        for (i = 0; i < lignes; i++) {
            for (j = 0; j < colonnes; j++) {
                a1[i * colonnes + j] = iDep.getPixelValue(j, i);
            }
        }

        for (i = 0; i < lignes; ++i) {
            icolonnes = i * colonnes;
            icol_1 = icolonnes - 1;
            icol_2 = icolonnes - 2;
            a2[icolonnes] = an1 * a1[icolonnes];
            a2[icolonnes + 1] = an1 * a1[icolonnes + 1]
                    + an2 * a1[icolonnes] - ad1 * a2[icolonnes];
            for (j = 2; j < colonnes; ++j) {
                a2[icolonnes + j] = an1 * a1[icolonnes + j] + an2 * a1[icol_1 + j]
                        - ad1 * a2[icol_1 + j] - ad2 * a2[icol_2 + j];
            }
        }

        for (i = 0; i < lignes; ++i) {
            icolonnes = i * colonnes;
            icol_1 = icolonnes + 1;
            icol_2 = icolonnes + 2;
            a3[icolonnes + col_1] = 0;
            a3[icolonnes + col_2] = an3 * a1[icolonnes + col_1];
            for (j = col_3; j >= 0; --j) {
                a3[icolonnes + j] = an3 * a1[icol_1 + j] + an4 * a1[icol_2 + j]
                        - ad1 * a3[icol_1 + j] - ad2 * a3[icol_2 + j];
            }
        }

        icol_1 = lignes * colonnes;

        for (i = 0; i < icol_1; ++i) {
            a2[i] += a3[i];
        }

        /*
         * FIRST STEP Y-GRADIENT : y-derivative
         */
        /*
         * columns top - downn
         */
        for (j = 0; j < colonnes; ++j) {
            a3[j] = 0;
            a3[colonnes + j] = an11 * a2[j] - ad1 * a3[j];
            for (i = 2; i < lignes; ++i) {
                a3[i * colonnes + j] = an11 * a2[(i - 1) * colonnes + j]
                        - ad1 * a3[(i - 1) * colonnes + j] - ad2 * a3[(i - 2) * colonnes + j];
            }
        }

        /*
         * columns down top
         */
        for (j = 0; j < colonnes; ++j) {
            a4[lig_1 * colonnes + j] = 0;
            a4[(lig_2 * colonnes) + j] = -an11 * a2[lig_1 * colonnes + j]
                    - ad1 * a4[lig_1 * colonnes + j];
            for (i = lig_3; i >= 0; --i) {
                a4[i * colonnes + j] = -an11 * a2[(i + 1) * colonnes + j]
                        - ad1 * a4[(i + 1) * colonnes + j] - ad2 * a4[(i + 2) * colonnes + j];
            }
        }

        icol_1 = colonnes * lignes;
        for (i = 0; i < icol_1; ++i) {
            a3[i] += a4[i];
        }

        for (i = 0; i < lignes; ++i) {
            for (j = 0; j < colonnes; ++j) {
                nf_gry[i * colonnes + j] = a3[i * colonnes + j];
            }
        }

        /*
         * SECOND STEP X-GRADIENT
         */
        for (i = 0; i < lignes; ++i) {
            for (j = 0; j < colonnes; ++j) {
                a1[i * colonnes + j] = (int) (iDep.getPixel(j, i));
            }
        }

        for (i = 0; i < lignes; ++i) {
            icolonnes = i * colonnes;
            icol_1 = icolonnes - 1;
            icol_2 = icolonnes - 2;
            a2[icolonnes] = 0;
            a2[icolonnes + 1] = an11 * a1[icolonnes];
            for (j = 2; j < colonnes; ++j) {
                a2[icolonnes + j] = an11 * a1[icol_1 + j]
                        - ad1 * a2[icol_1 + j] - ad2 * a2[icol_2 + j];
            }
        }

        for (i = 0; i < lignes; ++i) {
            icolonnes = i * colonnes;
            icol_1 = icolonnes + 1;
            icol_2 = icolonnes + 2;
            a3[icolonnes + col_1] = 0;
            a3[icolonnes + col_2] = -an11 * a1[icolonnes + col_1];
            for (j = col_3; j >= 0; --j) {
                a3[icolonnes + j] = -an11 * a1[icol_1 + j]
                        - ad1 * a3[icol_1 + j] - ad2 * a3[icol_2 + j];
            }
        }
        icol_1 = lignes * colonnes;
        for (i = 0; i < icol_1; ++i) {
            a2[i] += a3[i];
        }

        /*
         * on the columns
         */
        /*
         * columns top down
         */
        for (j = 0; j < colonnes; ++j) {
            a3[j] = an1 * a2[j];
            a3[colonnes + j] = an1 * a2[colonnes + j] + an2 * a2[j]
                    - ad1 * a3[j];
            for (i = 2; i < lignes; ++i) {
                a3[i * colonnes + j] = an1 * a2[i * colonnes + j] + an2 * a2[(i - 1) * colonnes + j]
                        - ad1 * a3[(i - 1) * colonnes + j] - ad2 * a3[(i - 2) * colonnes + j];
            }
        }

        /*
         * columns down top
         */
        for (j = 0; j < colonnes; ++j) {
            a4[lig_1 * colonnes + j] = 0;
            a4[lig_2 * colonnes + j] = an3 * a2[lig_1 * colonnes + j] - ad1 * a4[lig_1 * colonnes + j];
            for (i = lig_3; i >= 0; --i) {
                a4[i * colonnes + j] = an3 * a2[(i + 1) * colonnes + j] + an4 * a2[(i + 2) * colonnes + j]
                        - ad1 * a4[(i + 1) * colonnes + j] - ad2 * a4[(i + 2) * colonnes + j];
            }
        }

        icol_1 = colonnes * lignes;
        for (i = 0; i < icol_1; ++i) {
            a3[i] += a4[i];
        }

        for (i = 0; i < lignes; i++) {
            for (j = 0; j < colonnes; j++) {
                nf_grx[i * colonnes + j] = a3[i * colonnes + j];
            }
        }

        /*
         * SECOND STEP X-GRADIENT : the x-gradient is done
         */
        /*
         * THIRD STEP : NORM
         */
        /*
         * computatopn of the magnitude
         */
        for (i = 0; i < lignes; i++) {
            for (j = 0; j < colonnes; j++) {
                a2[i * colonnes + j] = nf_gry[i * colonnes + j];
            }
        }
        icol_1 = colonnes * lignes;
        for (i = 0; i < icol_1; ++i) {
            a2[i] = (float) Math.sqrt((a2[i] * a2[i]) + (a3[i] * a3[i]));
        }
        /*
         * THIRD STEP : the norm is done
         */
        result_array = new byte[nmem];

        //Recherche des niveaux min et max du gradiant
        double min = a2[0];
        double max = a2[0];
        for (i = 1; i < nmem; i++) {
            if (min > a2[i]) {
                min = a2[i];
            }
            if (max < a2[i]) {
                max = a2[i];
            }
        }

        //Normalisation de gradient de 0 a 255
        for (i = 0; i < nmem; ++i) {
            result_array[i] = (byte) (255 * (a2[i] / (max - min)));
        }

        //sauvegarde de la norme du gradiant
        iGrad.setPixels(result_array);

        return iGrad;
    }

    /**
     * main function for the snake
     *
     * @return Description of the Return Value
     */
    public double process() {
        int i;
        double force;
        Point2d displ = new Point2d();
        double maxforce = 0.0;
        double som = 0.0;
        double seuil = configuration.getGradThreshold();
        double DivForce = configuration.getMaxDisplacement();
        double minr = configuration.getRegMin();
        double maxr = configuration.getRegMax();
        double alpha = configuration.getAlpha();

        // EXPERIMENTAL
        double dist_plus = Prefs.get("ABSnake_ThreshDistPos.double", 100);
        double dist_minus = Prefs.get("ABSnake_ThreshDistNeg.double", 100);

        //IJ.log("process "+dist_plus+" "+dist_minus);

        for (i = 0; i < NPT; i++) {
            normale[i] = compute_normale(i);
        }
        block = 0;
        elimination = 0;
        for (i = 0; i < NPT; i++) {
            displ.x = 0.0;
            displ.y = 0.0;
            displ = searchTheClosestEdge(i, seuil, dist_plus, dist_minus, -1);

            force = Math.sqrt(Math.pow(displ.x, 2.0) + Math.pow(displ.y, 2.0));
            if (force > DivForce) {
                deplace[i].x = DivForce * (displ.x / force);
                deplace[i].y = DivForce * (displ.y / force);
            } else {
                deplace[i].x = displ.x;
                deplace[i].y = displ.y;
            }
            force = Math.sqrt(deplace[i].x * deplace[i].x + deplace[i].y * deplace[i].y);
            if (force > maxforce) {
                maxforce = force;
            }
            som += force;
        }
        dataDistance = som / NPT;

        for (i = 0; i < NPT; i++) {
            force = Math.sqrt(Math.pow(deplace[i].x, 2.0) + Math.pow(deplace[i].y, 2.0));
            lambda[i] = maxr / (1.0 + ((maxr - minr) / minr) * (force / maxforce));
        }
        if (elimination == 1) {
            destroysnake();
        }

        new_positions();
        resample(false);

        return dataDistance;
    }

    /**
     * SEGMENTATION : inside/outside the snake
     *
     * @param wi Description of the Parameter
     * @param he Description of the Parameter
     * @return binarised image (black=object inside snake)
     */
    public ByteProcessor segmentation(int wi, int he, int col) {
        Point2d pos;
        Point2d norm;
        Point2d ref;
        int top;
        int left;
        int right;
        int bottom;
        int i;
        int j;
        int x;
        int y;
        int count;
        double bden;
        double bnum;
        double bres;
        double lnorm;
        double ares;

        pos = new Point2d();
        norm = new Point2d();
        ref = new Point2d();

        ByteProcessor res = new ByteProcessor(wi, he);
        //res.invert();

        //Calcul des valeur permettant de definir le rectangle englobant du snake
        top = 0;
        bottom = 100000;
        left = 100000;
        right = 0;
        for (i = 0; i < NPT; i++) {
            if (points[i].y > top) {
                top = (int) points[i].y;
            }
            if (points[i].y < bottom) {
                bottom = (int) points[i].y;
            }
            if (points[i].x > right) {
                right = (int) points[i].x;
            }
            if (points[i].x < left) {
                left = (int) points[i].x;
            }
        }

        //On dessine l'interieur du snake a 255
        ref.x = 0;
        ref.y = 0;
        for (x = left; x < right; x++) {
            for (y = bottom; y < top; y++) {
                pos.x = x;
                pos.y = y;
                //norm.x = ref.x - pos.x;
                //norm.y = ref.y - pos.y;
                //lnorm = Math.sqrt(norm.x * norm.x + norm.y * norm.y);
                //norm.x /= lnorm;
                //norm.y /= lnorm;

                if (inside(pos)) {
                    res.putPixel(x, y, col);
                } else {
                    res.putPixel(x, y, 0);
                }
            }
        }
        return res;
    }

    /**
     * the point is inside the snake ?
     *
     * @param pos the point
     * @return inside ?
     */
    boolean inside(Point2d pos) {
        int count;
        int i;
        double bden;
        double bnum;
        double bres;
        double ares;
        double lnorm;
        Point2d norm = new Point2d();
        Point2d ref = new Point2d();

        ref.x = 0.0;
        ref.y = 0.0;
        norm.x = ref.x - pos.x;
        norm.y = ref.y - pos.y;
        lnorm = Math.sqrt(norm.x * norm.x + norm.y * norm.y);
        norm.x /= lnorm;
        norm.y /= lnorm;

        count = 0;
        for (i = 0; i < NPT - 1; i++) {
            bden = (-norm.x * points[i + 1].y + norm.x * points[i].y + norm.y * points[i + 1].x - norm.y * points[i].x);
            bnum = (-norm.x * pos.y + norm.x * points[i].y + norm.y * pos.x - norm.y * points[i].x);
            if (bden != 0) {
                bres = (bnum / bden);
            } else {
                bres = 5.0;
            }
            if ((bres >= 0.0) && (bres <= 1.0)) {
                ares = -(-points[i + 1].y * pos.x + points[i + 1].y * points[i].x
                        + points[i].y * pos.x + pos.y * points[i + 1].x - pos.y * points[i].x
                        - points[i].y * points[i + 1].x) / (-norm.x * points[i + 1].y
                        + norm.x * points[i].y + norm.y * points[i + 1].x - norm.y * points[i].x);
                if ((ares >= 0.0) && (ares <= lnorm)) {
                    count++;
                }
            }
        }
        // last point
        i = NPT - 1;
        bden = (-norm.x * points[0].y + norm.x * points[i].y + norm.y * points[0].x - norm.y * points[i].x);
        bnum = (-norm.x * pos.y + norm.x * points[i].y + norm.y * pos.x - norm.y * points[i].x);
        if (bden != 0) {
            bres = (bnum / bden);
        } else {
            bres = 5.0;
        }
        if ((bres >= 0.0) && (bres <= 1.0)) {
            ares = -(-points[0].y * pos.x + points[0].y * points[i].x
                    + points[i].y * pos.x + pos.y * points[0].x - pos.y * points[i].x
                    - points[i].y * points[0].x) / (-norm.x * points[0].y
                    + norm.x * points[i].y + norm.y * points[0].x - norm.y * points[i].x);
            if ((ares >= 0.0) && (ares <= lnorm)) {
                count++;
            }
        }
        return (count % 2 == 1);
    }
}
