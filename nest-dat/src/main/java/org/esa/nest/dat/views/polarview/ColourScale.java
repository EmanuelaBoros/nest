package org.esa.nest.dat.views.polarview;

import java.awt.*;
import java.awt.image.ColorModel;
import java.awt.image.IndexColorModel;
import java.util.Enumeration;
import java.util.Vector;

public class ColourScale {

    private ColorModel cm = null;
    private Color[] colors;
    private int thresholdCount;
    private int[] colorIndexThresholds;
    private final Vector<ColorBar> coloredClients = new Vector<ColorBar>();

    private double colorIndexValues[] = null;
    private double darkestValue = 0;
    private double darkestIndex = 0;
    private double lightestValue = 0;
    private double lightestIndex = 0;
    private double range = 0;
    private double scale = 0;

    ColourScale(Color colorTable[]) {
        thresholdCount = colorTable.length;
        colorIndexThresholds = new int[thresholdCount];
        colorIndexThresholds[0] = 0;
        colorIndexThresholds[thresholdCount - 1] = 255;
        colors = colorTable;
    }

    public static ColourScale newMonochromeScale(double range[], Color chromum) {
        final Color monochromeColorTable[] = {
                Color.black, Color.black, chromum, chromum
        };
        return new ColourScale(monochromeColorTable, range);
    }

    public static ColourScale newCustomScale(double range[]) {
        return new ColourScale(PolarView.waveColorTable, range);
    }

    ColourScale(Color colorTable[], double range[]) {
        this(colorTable);
        colorIndexValues = new double[thresholdCount];
        setRange(range[0], range[1]);
        createColorMap();
    }

    public void setRange(double range[]) {
        setRange(range[0], range[1]);
    }

    public void setRange(int range[]) {
        setRange(range[0], range[1]);
    }

    private void setRange(int minValue, int maxValue) {
        setRange((double) minValue, (double) maxValue);
    }

    private void setRange(double minValue, double maxValue) {
        darkestValue = minValue;
        lightestValue = maxValue;
        validateRange();
        setEvenThresholds();
        darkestIndex = colorIndexThresholds[0];
        lightestIndex = colorIndexThresholds[thresholdCount - 1];
        updateRange();
    }

    public byte getColorIndex(int value) {
        return getColorIndex((double) value);
    }

    public byte getColorIndex(float value) {
        return getColorIndex((double) value);
    }

    public byte getColorIndex(double value) {
        value -= darkestValue;
        if (value < 0.0D)
            return (byte) (int) darkestIndex;
        if (scale != 0.0D)
            value *= scale;
        value += darkestIndex;
        if (value > lightestIndex)
            return (byte) (int) lightestIndex;
        else
            return (byte) ((int) Math.round(value) & 0xff);
    }

    private int getIntegerColorValue(int index) {
        return (int) Math.round(getDoubleColorValue(index));
    }

    private float getFloatColorValue(int index) {
        return (float) getDoubleColorValue(index);
    }

    private double getDoubleColorValue(int index) {
        double value = (double) index - darkestIndex;
        if (scale != 0.0D)
            value /= scale;
        return value + darkestValue;
    }

    private int getIntegerThresholdValue(int thresholdIndex) {
        return (int) Math.round(getDoubleThresholdValue(thresholdIndex));
    }

    private float getFloatThresholdValue(int thresholdIndex) {
        return (float) getDoubleThresholdValue(thresholdIndex);
    }

    private double getDoubleThresholdValue(int thresholdIndex) {
        return colorIndexValues[thresholdIndex];
    }

    private void updateColorValues() {
        for (int i = 0; i < thresholdCount; i++) {
            colorIndexValues[i] = getIntegerColorValue(colorIndexThresholds[i]);
        }
    }

    private void validateRange() {
        darkestValue = Math.min(darkestValue, lightestValue);
        range = lightestValue - darkestValue;
        scale = 255D / range;
    }

    public final Color getColor(int value) {
        return new Color(getRGB(value));
    }

    public final Color getColor(float value) {
        return new Color(getRGB(value));
    }

    public final Color getColor(double value) {
        return new Color(getRGB(value));
    }

    private int getRGB(int value) {
        return cm.getRGB(getColorIndex(value) & 0xff);
    }

    private int getRGB(float value) {
        return cm.getRGB(getColorIndex(value) & 0xff);
    }

    private int getRGB(double value) {
        return cm.getRGB(getColorIndex(value) & 0xff);
    }

    public ColorModel getColorModel() {
        return cm;
    }

    public synchronized void addColoredObject(ColorBar ip) {
        if (!coloredClients.contains(ip)) {
            coloredClients.addElement(ip);
        }
    }

    private int limitColorThreshold(int i, int index) {
        final int lastThreshold = thresholdCount - 1;
        if (i < 1 || i >= lastThreshold)
            return index;
        int upperLimit = colorIndexThresholds[i + 1];
        if (i < lastThreshold - 1)
            upperLimit -= 2;
        int lowerLimit = colorIndexThresholds[i - 1];
        if (i > 1)
            lowerLimit += 2;
        int idx = Math.min(index, upperLimit);
        idx = Math.max(idx, lowerLimit);
        return idx;
    }

    private void setEvenThresholds() {
        int N = thresholdCount - 1;
        int first = 0;
        int last = N;
        if (colors[last].equals(colors[last - 1])) {
            colorIndexThresholds[last] = 255;
            last--;
            N--;
        }
        if (colors[first].equals(colors[first + 1])) {
            colorIndexThresholds[first] = 0;
            first++;
            N--;
        }
        final double colorStep = 255D / N;
        int offset = 0;
        int i = 0;
        for (int t = first; t <= last; t++) {
            colorIndexThresholds[t] = offset + (int) Math.round((double) i * colorStep);
            i++;
        }
    }

    private void createColorMap() {
        final int lastThreshold = thresholdCount - 1;
        final byte cmap[] = new byte[768];
        Color lastColor = colors[0];
        int lastIndex = colorIndexThresholds[0];
        int c = 0;
        for (int i = 1; i < thresholdCount; i++) {
            final int cRange = colorIndexThresholds[i] - lastIndex;
            final int lastRGB[] = {lastColor.getRed(), lastColor.getGreen(), lastColor.getBlue()};
            final int nextRGB[] = {colors[i].getRed(), colors[i].getGreen(), colors[i].getBlue()};

            for (int j = 0; j < cRange; j++) {
                final float nextScale = (float) j / (float) cRange;
                final float lastScale = 1.0F - nextScale;
                cmap[c++] = (byte) (int) ((float) lastRGB[0] * lastScale + (float) nextRGB[0] * nextScale);
                cmap[c++] = (byte) (int) ((float) lastRGB[1] * lastScale + (float) nextRGB[1] * nextScale);
                cmap[c++] = (byte) (int) ((float) lastRGB[2] * lastScale + (float) nextRGB[2] * nextScale);
            }

            lastColor = colors[i];
            lastIndex = colorIndexThresholds[i];
        }

        final Color finalColor = colors[lastThreshold];
        cmap[c++] = (byte) finalColor.getRed();
        cmap[c++] = (byte) finalColor.getGreen();
        cmap[c] = (byte) finalColor.getBlue();
        cm = new IndexColorModel(8, 256, cmap, 0, false);
    }

    private synchronized void notifyMapChange() {
        ColorBar ip;
        for (Enumeration elem = coloredClients.elements(); elem.hasMoreElements(); ip.updatedColorMap()) {
            ip = (ColorBar) elem.nextElement();
        }
    }

    private synchronized void notifyRangeChange() {
        ColorBar ip;
        for (Enumeration elem = coloredClients.elements(); elem.hasMoreElements(); ip.updatedColorScale()) {
            ip = (ColorBar) elem.nextElement();
        }
    }

    private void updateRange() {
        updateColorValues();
        createColorMap();
        notifyRangeChange();
    }

    private void update() {
        updateColorValues();
        createColorMap();
        notifyMapChange();
    }
}