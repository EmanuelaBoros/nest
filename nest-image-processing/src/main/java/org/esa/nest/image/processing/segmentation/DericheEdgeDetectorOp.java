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

import com.bc.ceres.core.ProgressMonitor;
import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.MultiLineString;
import com.vividsolutions.jts.geom.Polygon;
import com.vividsolutions.jts.geom.PrecisionModel;
import com.vividsolutions.jts.util.GeometricShapeFactory;
import ij.ImagePlus;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;
import java.awt.Rectangle;
import java.awt.image.BufferedImage;
import java.awt.image.RenderedImage;
import java.util.HashMap;
import java.util.Map;
import org.esa.beam.framework.datamodel.Band;
import org.esa.beam.framework.datamodel.Product;
import org.esa.beam.framework.datamodel.ProductData;
import org.esa.beam.framework.gpf.OperatorException;
import org.esa.beam.framework.gpf.OperatorSpi;
import org.esa.beam.framework.gpf.Tile;
import org.esa.beam.framework.gpf.annotations.OperatorMetadata;
import org.esa.beam.framework.gpf.annotations.Parameter;
import org.esa.beam.framework.gpf.annotations.SourceProduct;
import org.esa.beam.framework.gpf.annotations.TargetProduct;
import org.esa.nest.gpf.OperatorUtils;
import ij.measure.*;
import ij.gui.*;
import ij.*;
import ij.plugin.frame.PlugInFrame;
import ij.process.ImageStatistics;
import java.awt.*;
import java.awt.event.*;
import java.awt.image.ImageProducer;
import java.net.URL;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import javax.swing.JOptionPane;
import org.esa.beam.framework.datamodel.PlainFeatureFactory;
import org.esa.beam.framework.datamodel.ProductNodeGroup;
import org.esa.beam.framework.datamodel.VectorDataNode;
import org.esa.beam.framework.gpf.Operator;
import org.esa.beam.jai.ImageManager;
import org.geotools.feature.FeatureCollection;
import org.geotools.feature.FeatureIterator;
import org.geotools.geometry.GeometryBuilder;
import org.opengis.feature.simple.SimpleFeature;
import org.opengis.feature.simple.SimpleFeatureType;
import org.opengis.referencing.crs.CoordinateReferenceSystem;
import org.openimaj.image.FImage;
import org.openimaj.image.ImageUtilities;
import org.openimaj.image.processing.edges.CannyEdgeDetector;

/**
 * This plug-in takes as parameters a grayscale image and two thresholds (low
 * and high), and returns the DericheEdgeDetectorOp image
 *
 * @author Emanuela Boros
 */
@OperatorMetadata(alias = "DericheEdgeDetector",
category = "SAR Tools\\Image Processing",
description = "DericheEdgeDetector")
public class DericheEdgeDetectorOp extends Operator {

    public static float[] probabilityHistogram;
    final static int MAX_VALUE = 256;
    final static int MIN_VALUE = 0;
    public static int N;
    @SourceProduct(alias = "source")
    private Product sourceProduct = null;
    @TargetProduct
    private Product targetProduct;
    @Parameter(description = "The list of source bands.", alias = "sourceBands", itemAlias = "band",
    rasterDataNodeType = Band.class, label = "Source Bands")
    private String[] sourceBandNames;
    @Parameter(description = "HighThreshold", defaultValue = "100", label = "HighThreshold")
    private float highThreshold = 100f;
    @Parameter(description = "LowThreshold", defaultValue = "10", label = "LowThreshold")
    private float lowThreshold = 10f;
    private final Map<String, String[]> targetBandNameToSourceBandName = new HashMap<String, String[]>();
    private int sourceImageWidth;
    private int sourceImageHeight;
    private static boolean processed = false;
    private int halfSizeX;
    private int halfSizeY;
    private int filterSizeX = 3;
    private int filterSizeY = 3;
    private static ImagePlus fullImagePlus;
//    private static ByteProcessor fullByteProcessor;
    private static CannyEdgeDetector cannyEdgeDetector;
    /**
     *
     *
     */
    private ImagePlus img, prebleach, postbleach; // holds the prebleach-postbleach slice
    private ImageStack stack;
    private ImageWindow current;
    private ImageProcessor ippre, ippost; // image processors for stack and slices
    private int preref = 1; // holds the slice number for the prebleach slice
    private int postref = 2; // holds the slice number for the postbleach slice
    private Region frap, ref, base, whole;
    private float alphaD = (float) 0.5;
    private float upperTreshold = (float) 50;
    private float lowerTreshold = (float) 100;
    private double radius = 5;
    private int nSlices = 1;
    private ResultsTable rtable = new ResultsTable(), normtable = new ResultsTable();
    private boolean doublenorm = true, autocalc = false, image = false;
    private int[] time;
    private int interval = 1;
    float wholepre, frappre;

    /**
     * Initializes this operator and sets the one and only target product.
     * <p>The target product can be either defined by a field of type
     * {@link org.esa.beam.framework.datamodel.Product} annotated with the
     * {@link org.esa.beam.framework.gpf.annotations.TargetProduct TargetProduct}
     * annotation or by calling {@link #setTargetProduct} method.</p> <p>The
     * framework calls this method after it has created this operator. Any
     * client code that must be performed before computation of tile data should
     * be placed here.</p>
     *
     * @throws org.esa.beam.framework.gpf.OperatorException If an error occurs
     * during operator initialization.
     * @see #getTargetProduct()
     */
    @Override
    public void initialize() throws OperatorException {

        cannyEdgeDetector = new CannyEdgeDetector();

        try {
            sourceImageWidth = sourceProduct.getSceneRasterWidth();
            sourceImageHeight = sourceProduct.getSceneRasterHeight();

            halfSizeX = filterSizeX / 2;
            halfSizeY = filterSizeY / 2;

            createTargetProduct();

        } catch (Throwable e) {
            OperatorUtils.catchOperatorException(getId(), e);
        }
    }

    /**
     * Create target product.
     *
     * @throws Exception The exception.
     */
    private void createTargetProduct() throws Exception {

        targetProduct = new Product(sourceProduct.getName(),
                sourceProduct.getProductType(),
                sourceImageWidth,
                sourceImageHeight);

        OperatorUtils.copyProductNodes(sourceProduct, targetProduct);

        OperatorUtils.addSelectedBands(
                sourceProduct, sourceBandNames, targetProduct, targetBandNameToSourceBandName, true);
    }

    /**
     * Called by the framework in order to compute a tile for the given target
     * band. <p>The default implementation throws a runtime exception with the
     * message "not implemented".</p>
     *
     * @param targetBand The target band.
     * @param targetTile The current tile associated with the target band to be
     * computed.
     * @param pm A progress monitor which should be used to determine
     * computation cancelation requests.
     * @throws org.esa.beam.framework.gpf.OperatorException If an error occurs
     * during computation of the target raster.
     */
    @Override
    public void computeTile(Band targetBand, Tile targetTile, ProgressMonitor pm) throws OperatorException {

        try {
            final Rectangle targetTileRectangle = targetTile.getRectangle();
            final int x0 = targetTileRectangle.x;
            final int y0 = targetTileRectangle.y;
            final int w = targetTileRectangle.width;
            final int h = targetTileRectangle.height;

            final Rectangle sourceTileRectangle = getSourceTileRectangle(x0, y0, w, h);
            Tile sourceRaster;
            final String[] srcBandNames = targetBandNameToSourceBandName.get(targetBand.getName());
            Band sourceBand = sourceProduct.getBand(srcBandNames[0]);
            sourceRaster = getSourceTile(sourceBand, sourceTileRectangle);
            if (sourceRaster == null) {
                throw new OperatorException("Cannot get source tile");
            }

            computeDericheEdgeDetection(sourceBand, sourceRaster, targetTile, x0, y0, w, h, pm);

        } catch (Throwable e) {
            OperatorUtils.catchOperatorException(getId(), e);
        } finally {
            pm.done();
        }
    }
    /**
     * Apply HysteresisThesholding
     *
     * @param sourceRaster The source tile for the band.
     * @param targetTile The current tile associated with the target band to be
     * computed.
     * @param x0 X coordinate for the upper-left point of the
     * target_Tile_Rectangle.
     * @param y0 Y coordinate for the upper-left point of the
     * target_Tile_Rectangle.
     * @param w Width for the target_Tile_Rectangle.
     * @param h Hight for the target_Tile_Rectangle.
     * @param pm A progress monitor which should be used to determine
     * computation cancellation requests.
     * @throws org.esa.beam.framework.gpf.OperatorException If an error occurs
     * during computation of the filtered value.
     */
    private static FImage fImage;
    private static BufferedImage bf;
    private static BufferedImage fullBufferedImage;
    private static ArrayList<Roi> currentROIs = new ArrayList<Roi>();

    private synchronized void computeDericheEdgeDetection(final Band sourceBand, final Tile sourceRaster,
            final Tile targetTile, final int x0, final int y0, final int w, final int h,
            final ProgressMonitor pm) {

        if (!isProcessed()) {
            final RenderedImage fullRenderedImage = sourceBand.getSourceImage().getImage(0);
            fullBufferedImage = new BufferedImage(sourceBand.getSceneRasterWidth(),
                    sourceBand.getSceneRasterHeight(),
                    BufferedImage.TYPE_USHORT_GRAY);
            fullBufferedImage.setData(fullRenderedImage.getData());
            fullImagePlus = new ImagePlus(sourceBand.getDisplayName(), fullBufferedImage);

            fImage = ImageUtilities.createFImage(fullBufferedImage);
//            cannyEdgeDetector.processImage(fImage);

            ProductNodeGroup<VectorDataNode> productNodeGroup = sourceProduct.getVectorDataGroup();
            for (int i = 0; i < productNodeGroup.getNodeCount(); i++) {
                FeatureCollection<SimpleFeatureType, SimpleFeature> features =
                        productNodeGroup.get(i).getFeatureCollection();
                FeatureIterator<SimpleFeature> featureIterator = features.features();

                while (featureIterator.hasNext()) {
                    SimpleFeature feature = featureIterator.next();
                    Object value = feature.getDefaultGeometry();
                    Polygon polygon = (Polygon) value;
                    int xx[] = new int[polygon.getCoordinates().length];
                    int yy[] = new int[polygon.getCoordinates().length];

                    for (int j = 0; j < polygon.getCoordinates().length; j++) {
                        Coordinate coordinate = polygon.getCoordinates()[j];
                        xx[j] = (int) (coordinate.x);
                        yy[j] = (int) (coordinate.y);
                    }

                    PolygonRoi currentROI = new PolygonRoi(xx, yy,
                            polygon.getCoordinates().length - 1, Roi.FREEROI);
                    currentROIs.add(new Roi(currentROI.getBounds()));

                }
//                JOptionPane.showMessageDialog(null, polygons.size() + " Rois",
//                        "getImagePlus", JOptionPane.INFORMATION_MESSAGE);
            }

            if (currentROIs.size() > 0) {
                final ImageProcessor fullImageProcessor = fullImagePlus.getProcessor();

                ByteProcessor fullByteProcessor = (ByteProcessor) fullImageProcessor.convertToByte(true);

                for (int i = 0; i < currentROIs.size(); i++) {

                    Roi currentROI = currentROIs.get(i);
                    fullByteProcessor.setRoi(new Roi(currentROI.getBounds()));

                    ImageProcessor roiImageProcessor = fullByteProcessor.crop();

                    ImagePlus currentImagePlus = new ImagePlus(sourceBand.getName() + "#" + i);
                    currentImagePlus.setProcessor(roiImageProcessor);
                    currentImagePlus.show();
                }
            }
//            saveToGeometry(polys, sourceBand.getName());

            bf = new BufferedImage(sourceBand.getSceneRasterWidth(),
                    sourceBand.getSceneRasterHeight(),
                    BufferedImage.TYPE_USHORT_GRAY);

//            bf = ImageUtilities.createWorkingImage(ImageUtilities.createBufferedImage(fImage));
            bf = ImageUtilities.createBufferedImageForDisplay(fImage);

//            try {
//                ImageUtilities.write(fImage, new File("D:/final.png"));
//            } catch (IOException ex) {
//                Logger.getLogger(DericheEdgeDetectorOp.class.getName()).log(Level.SEVERE, null, ex);
//            }

//            fullImagePlus = new ImagePlus(sourceBand.getDisplayName(), fullBufferedImage);
//
//            final ImageProcessor fullImageProcessor = fullImagePlus.getProcessor();
//
//            fullByteProcessor = (ByteProcessor) fullImageProcessor.convertToByte(true);

//            ImagePlus imagePlus = new ImagePlus(sourceBand.getDisplayName()
//                    + y0 + x0 + w + h, fullByteProcessor);
//            imagePlus.setProcessor(fullByteProcessor);

//            imagePlus = new Image_Edge().run(imagePlus);

//            fullByteProcessor = (ByteProcessor) imagePlus.getProcessor().convertToByte(true);
            setProcessed(true);
        }

        final Rectangle srcTileRectangle = sourceRaster.getRectangle();
        BufferedImage newBf = bf.getSubimage(srcTileRectangle.x, srcTileRectangle.y,
                srcTileRectangle.width, srcTileRectangle.height);
        FImage crop = ImageUtilities.createFImage(newBf);

//                fImage.extractROI(srcTileRectangle.x, srcTileRectangle.y,
//                srcTileRectangle.width, srcTileRectangle.height);
////        if (null != fImage) {
//            crop = fImage.extractROI(srcTileRectangle.x, srcTileRectangle.y,
//                    srcTileRectangle.width, srcTileRectangle.height);
////        }

//        ImageProcessor aPartProcessor = fullByteProcessor.duplicate();

//        aPartProcessor.setRoi(srcTileRectangle);


//        ImageProcessor roiImageProcessor = aPartProcessor.crop();
//        roiImageProcessor = (ByteProcessor) roiImageProcessor.convertToByte(true);

        final ProductData trgData = targetTile.getDataBuffer();
        final ProductData sourceData =
                ProductData.createInstance(crop.toPackedARGBPixels());
//        ProductData.createInstance((byte[]) roiImageProcessor.getPixels());

        final int maxY = y0 + h;
        final int maxX = x0 + w;
        for (int y = y0; y < maxY; ++y) {
            for (int x = x0; x < maxX; ++x) {

                trgData.setElemFloatAt(targetTile.getDataBufferIndex(x, y),
                        sourceData.getElemFloatAt(sourceRaster.getDataBufferIndex(x, y)));
            }
        }
    }

//    private static Rectangle createRectangle(Coordinate centre,
//            double width, double height, PrecisionModel precision) {
//
//        GeometricShapeFactory gsf = new GeometricShapeFactory(
//                getGeometryFactory(precision));
//
//        gsf.setCentre(centre);
//        gsf.setWidth(width);
//        gsf.setHeight(height);
//        gsf.setNumPoints(4);
//        return gsf.createRectangle();
//    }
    public Rectangle getRectangle(Polygon polygon) {

        Rectangle rect = new Rectangle(Integer.MAX_VALUE, Integer.MIN_VALUE, 0, 0);
        for (Coordinate p : polygon.getCoordinates()) {
            rect.x = (int) Math.min(p.x, rect.x);
            rect.y = (int) Math.min(p.y, rect.y);
        }
//        for (Coordinate p : polygon.getCoordinates()) {
//            rect.right((int) Math.max(p.x, rect.right()));
//            rect.bottom((int) Math.max(p.y, rect.bottom()));
//        }
        return rect;
    }
    /**
     * A static reference to the JTS-Double-Precision-Model
     */
    public static PrecisionModel DOUBLE_PRECISIONMODEL =
            new PrecisionModel(PrecisionModel.FLOATING);
    /**
     * A static reference to the JTS-Float-Precision-Model
     */
    public static PrecisionModel FLOAT_PRECISIONMODEL =
            new PrecisionModel(PrecisionModel.FLOATING_SINGLE);
    private static GeometryFactory DOUBLE_GEOMETRYFACTORY =
            new GeometryFactory(DOUBLE_PRECISIONMODEL);
    private static GeometryFactory FLOAT_GEOMETRYFACTORY =
            new GeometryFactory(FLOAT_PRECISIONMODEL);

    private static GeometryFactory getGeometryFactory(PrecisionModel p) {
        return isDouble(p) ? DOUBLE_GEOMETRYFACTORY
                : isFloat(p) ? FLOAT_GEOMETRYFACTORY
                : new GeometryFactory(p);
    }

    private static boolean isFloat(PrecisionModel p) {
        return p.isFloating() && p.getMaximumSignificantDigits() == 6;
    }

    private static boolean isDouble(PrecisionModel p) {
        return p.isFloating() && p.getMaximumSignificantDigits() == 16;
    }

    private void saveToGeometry(final Collection<MultiLineString> polys, final String bandName) {

        final CoordinateReferenceSystem modelCrs = ImageManager.getModelCrs(
                targetProduct.getGeoCoding());
        final SimpleFeatureType type = PlainFeatureFactory.createDefaultFeatureType(modelCrs);

        final VectorDataNode vectorDataNode = new VectorDataNode(bandName + "_vector", type);
        final FeatureCollection<SimpleFeatureType, SimpleFeature> collection =
                vectorDataNode.getFeatureCollection();

        final String style = vectorDataNode.getDefaultStyleCss();

        final Iterator<MultiLineString> iter = polys.iterator();
        int k = 0;
        while (k < polys.size()) {
            final MultiLineString p = iter.next();
            p.normalize();

            final String name = bandName + '_' + p.getUserData() + '_' + k;

            final SimpleFeature feature = PlainFeatureFactory.createPlainFeature(type, name, p, style);

            collection.add(feature);

            // print the polygon as a WKT string
            //System.out.println(p.toText());
            k++;
        }

        targetProduct.getVectorDataGroup().add(vectorDataNode);
    }

    public void run(ImagePlus image) {
//        if (IJ.versionLessThan("1.36")) {
//            IJ.error("Version 1.36 or higher required.\n"
//                    + "Please update your ImageJ."); // make sure right version is used
//        }
//        if (arg.equals("about")) {
//            showAbout();
//            return;
//        } // end if
//
//        if (arg.equals("Pic")) {
//        showPic(image);
//            return;
//        } else {
//            image = getNewImage();

        nSlices = image.getStackSize();
        if (nSlices < 2) {
            IJ.error("Stack required");
        } // end if
        time = new int[nSlices];

        frap = new Region("FRAP region", nSlices);
        ref = new Region("Reference", nSlices);
        base = new Region("Background", nSlices);
        whole = new Region("Whole cell", nSlices);
        Display display = new Display(); // creating new frame
        GUI.center(display);
        display.setVisible(true);

//            if (image) {
        displayGuidance(image);
//            }
//        }

    } // END run FRAP_Analysis

    /**
     * Searches for the current open image. Returns true if found false if
     * otherwise
     *
     */
    private boolean getNewImage() {
        img = WindowManager.getCurrentImage();
        current = WindowManager.getCurrentWindow();
        if (img == null) {
            IJ.showMessage("No images open.");
            return false;
        }
        nSlices = img.getStackSize();
        if (nSlices < 2) {
            IJ.error("Stack required");
            return false;
        } // end if
        time = new int[nSlices];
        return true;

    }

    /**
     * Displays the guidance pictures for fast identification of the different
     * regions
     */
    private void displayGuidance(ImagePlus img) {
//        stack = img.getStack();
        ippre = ImageEdge.areaEdge(img.getProcessor(),
                radius, alphaD, upperTreshold, lowerTreshold);
        prebleach = new ImagePlus("Pre Bleach", ippre);
        prebleach.show();

        ippost = ImageEdge.areaEdge(img.getProcessor(),
                radius, alphaD, upperTreshold, lowerTreshold);
        postbleach = new ImagePlus("Post Bleach", ippost);
        postbleach.show();
    } // END displayGuidance

    /**
     * @return the processed
     */
    public synchronized boolean isProcessed() {
        return processed;
    }

    /**
     * @param processed the processed to set
     */
    public synchronized void setProcessed(boolean processed) {
        this.processed = processed;
    }

    @SuppressWarnings("serial")
    /**
     * Creates the interface
     *
     */
    private class Display extends PlugInFrame implements ActionListener, ItemListener {

        Frame window;
        Panel normpanel, autopanel, resetpanel; // for the normalization menu
        Choice norm, auto;

        public Display() {

            super("FRAP Analysis"); // call construction for the PlugInFrame
            if (window != null) {	// prevents construction twice
                window.toFront();
                return;
            } // end if

            window = this;
            setLayout(new GridLayout(10, 1, 5, 5));

            // choice menus

            normpanel = new Panel();
            Choice norm = new Choice();
            normpanel.add(new Label("Normalization :"));
            norm.add("Double");
            norm.add("Single");
            norm.addItemListener(this);
            normpanel.add(norm);

            autopanel = new Panel();
            Choice auto = new Choice();
            autopanel.add(new Label("Measurement :"));
            auto.add("Manual");
            auto.add("Automatic");
            auto.addItemListener(this);
            autopanel.add(auto);

            // reset panel
            resetpanel = new Panel();
            resetpanel.setLayout(new GridLayout(1, 2));
            Button resetbutton = new Button("Reset");
            resetbutton.addActionListener(this);
            Button imagebutton = new Button("New image");
            imagebutton.addActionListener(this);
            resetpanel.add(resetbutton);
            resetpanel.add(imagebutton);

            // window construction
            window.add(frap.panel);
            window.add(base.panel);
            window.add(whole.panel);
            window.add(ref.panel);
            add(autopanel);
            add(normpanel);
            addButton("Normalize");
            add(resetpanel);
            addButton("Settings");
            addButton("Help");
            pack();

        } // END constructor Display

        public void windowClosing(WindowEvent e) {
            super.windowClosing(e);
            window = null; // resets the window
        } // END windowClosing Display

        private void addButton(String label) {
            Button b = new Button(label);
            b.addActionListener(this);
            window.add(b);
        } // END addButton Display

        public void actionPerformed(ActionEvent e) {
            String action = e.getActionCommand();
            if (action.equals("Normalize")) {
                if (image) {
                    Normalize();
                }
            } else if (action.equals("Reset")) {
                Reset();
            } else if (action.equals("Settings")) {
                Settings();
            } else if (action.equals("Help")) {
                showAbout();
            } else if (action.equals("New image")) {
                if (image) {
                    prebleach.close();
                    postbleach.close();
                }
                image = getNewImage();
                if (image) {
                    Reset();
                    displayGuidance(img);
                }

            }
        } // END actionPerformed Display

        public void itemStateChanged(ItemEvent e) {
            String item = (String) e.getItem();
            if (item.equals("Double")) {
                doublenorm = true;
            }
            if (item.equals("Single")) {
                doublenorm = false;
            }
            if (item.equals("Manual")) {
                autocalc = false;
            }
            if (item.equals("Automatic")) {
                autocalc = true;
            }
        }
    } // END Display

    /**
     * Sets up the arrays and the methods for the different regions, and links
     * the buttons on the panel to those regions.
     *
     */
    private class Region implements ActionListener {

        Panel panel;
        Button measurebutton, applybutton, setbutton;
        Roi roi;
        float[] intensity, area;
        int[] slice;
        String label;

        public Region(String label, int length) {
            this.label = label;
            this.intensity = new float[length];
            this.area = new float[length];
            this.slice = new int[length];
            panel = new Panel();
            panel.setLayout(new GridLayout(1, 4));
            panel.add(new Label(label + " :"));

            setbutton = new Button("Set ROI");
            setbutton.addActionListener(this);
            applybutton = new Button("Apply to Image");
            applybutton.addActionListener(this);
            measurebutton = new Button("Measure");
            measurebutton.addActionListener(this);
            panel.add(setbutton);
            panel.add(applybutton);
            panel.add(measurebutton);
            roi = null;
        } // END constructor Region

        public void actionPerformed(ActionEvent e) {
            String action = e.getActionCommand();
            if (image) {
                if (action.equals("Set ROI")) {
                    this.roi = this.roiSet();
                } else if (action.equals("Apply to Image")) {
                    this.Apply(roi);
                } else if (action.equals("Measure")) {
                    this.imgMeasure();
                }
            } else {
                IJ.showMessage("No images open. \nOpen an image and click \"New image\".");
            }
        } // END actionPerformed Region

        private Roi roiSet() {
            ImagePlus img = getImage();
            if (img == null) {
                return null;
            }
            Roi roi = img.getRoi();
            if (roi == null) {
                IJ.showMessage("No ROI selected");
                return null;
            }
            return roi;
        }

        private void Apply(Roi roi) {
            WindowManager.setCurrentWindow(current);
            ImagePlus img = getImage();
            if (img == null) {
                return;
            } else if (roi == null) {
                IJ.showMessage("No ROI available. Set ROI first");
                return;
            }
            img.setRoi(roi);
        }

        private void imgMeasure() {
            ImagePlus img = getImage();
            if (img == null) {
                return;
            }
            Roi roi = img.getRoi();
            if (autocalc) {
                for (int i = 1; i <= nSlices; i++) {
                    img.setSlice(i);
                    Measure(img, roi);
                }
                img.setSlice(1);
                return;
            } else {
                Measure(img, roi);
            }
        }

        private void Measure(ImagePlus slice, Roi roi) {
            ImageProcessor ip = img.getProcessor();
            ip.setRoi(roi);
            ImageStatistics stats = ImageStatistics.getStatistics(ip, Measurements.MEAN + Measurements.AREA, null);
            int index = img.getCurrentSlice() - 1;
            this.slice[index] = index + 1;
            this.area[index] = (float) stats.area;
            this.intensity[index] = (float) stats.mean;

            displayResults(label, index + 1, stats.area, stats.mean);
        }

        private ImagePlus getImage() {
            ImagePlus imp = WindowManager.getCurrentImage();
            if (imp == null) {
                IJ.showMessage("There are no images open.");
                return null;
            } else {
                return imp;
            }
        }
    } // END private abstract class Region

    /**
     * Allows to change the settings of the plugin
     */
    private void Settings() {
        GenericDialog gd = new GenericDialog("Parameters");
        gd.addNumericField("median filter radius", radius, 2);
        gd.addNumericField("Deriche alpha value", alphaD, 2);
        gd.addNumericField("Hysteresis High threshold", upperTreshold, 2);
        gd.addNumericField("Hysteresis Low threshold", lowerTreshold, 2);
        gd.addNumericField("Pre Bleaching slice", preref, 0);
        gd.addNumericField("Post Bleaching slice", postref, 0);
        gd.addNumericField("Time units between slides", interval, 0);
        gd.showDialog();
        radius = gd.getNextNumber();
        alphaD = (float) gd.getNextNumber();
        upperTreshold = (float) gd.getNextNumber();
        lowerTreshold = (float) gd.getNextNumber();
        preref = (int) gd.getNextNumber();
        postref = (int) gd.getNextNumber();
        interval = (int) gd.getNextNumber();
        Reset();
    }

    public void displayResults(String label, int slice, double area, double mean) {
        rtable.incrementCounter();
        rtable.addLabel("Region", label);
        rtable.addValue("Slice", slice);
        rtable.addValue("Area", area);
        rtable.addValue("Intensity", mean);
        rtable.show("Measurements");
    }

    public void Normalize() {

        float[] frapnorm = new float[nSlices];
        frappre = 0;
        wholepre = 0;
        normtable.reset();
        // set time intervals
        for (int i = 0; i < nSlices; i++) {
            time[i] = (i + 1 - postref) * interval;
        }

        // Calculate average intensities prebleach
        for (int i = 0; i < postref - 1; i++) {
            wholepre += (whole.intensity[i] - base.intensity[i]);
            frappre += (frap.intensity[i] - base.intensity[i]);
        }
        wholepre = wholepre / preref;
        frappre = frappre / preref;

        // Calculate single normalization
        for (int i = 0; i < nSlices; i++) {
            frapnorm[i] = (frap.intensity[i] - base.intensity[i]) / frappre;
        }
        if (doublenorm) {
            for (int i = 0; i < nSlices; i++) {
                frapnorm[i] = frapnorm[i] * (wholepre / (whole.intensity[i] - base.intensity[i]));
            }
        }

        // display results
        for (int i = 0; i < nSlices; i++) {
            normtable.incrementCounter();
            normtable.addLabel(" ", "Intensities for slice " + (i + 1) + " :");
            normtable.addValue("Normalized", frapnorm[i]);
            normtable.addValue("Frap", frap.intensity[i]);
            normtable.addValue("Whole", whole.intensity[i]);
            normtable.addValue("Base", base.intensity[i]);
            normtable.addValue("time", time[i]);
            if (ref.intensity != null) {
                normtable.addValue("Ref", ref.intensity[i]);
            }

        }
        normtable.show("Normalization Results");
    }

    public void Reset() {
        setNull(frap);
        setNull(base);
        setNull(ref);
        setNull(whole);
        rtable.reset();
    }

    private void setNull(Region reg) {
        reg.area = new float[nSlices];
        reg.intensity = new float[nSlices];
        reg.slice = new int[nSlices];


    }

    private void showAbout() { // About message
        IJ.showMessage("About FRAP_Norm",
                "This plugin is specifically constructed to extract data from stacks of images\n"
                + "for FRAP analysis. It implements a median filter, deriche filter and subsequent\n"
                + "hysteresis to find the outline of the different regions in the cell. Parameters\n"
                + "for this filters can be put in the settings, together with the time interval and\n"
                + "the slice numbers for the prebleach and postbleach guidance images.\n \n"
                + "It is written originally for analysis of chromatin dynamics in  Arabidopsis thaliana,\n"
                + "and uses the normalization methods outlined in Phair et al (2004):\n"
                + "Measurement of dynamic protein binding to chromatin in vivo using photobleaching \n"
                + "microscopy, Methods Enzymol 375, 393-414.\n"
                + "Both single and double normalization can be carried out.\n \n"
                + "Procedure :\n"
                + "Use the ROI-tools, e.g. the wand tool or the rectangular tool to select\n"
                + "a region of interest. Use the \"Set ROI\" button next to the corresponding region\n"
                + "to save the ROI for that region. After defining the different ROIs, select the slice\n"
                + "you want to measure and press \"Apply to Image\" for a region. This will show up\n"
                + "the ROI for that region. Adjust the ROI by dragging it if necessary, and press\n"
                + "the \"Measure\" button. This will show up the measurement in the results window.\n \n"
                + "When all measurements are done, you can press the \"Normalize\" button to normalize\n"
                + "the measurements. Normalization is done as outlined in Phair et al(2004), and can be\n"
                + "done with (double) or without (single) the use of the whole cell measurements.\n"
                + "Measurement of the reference area is not obliged for normalization.\n \n"
                + "Use the \"Reset\" button to erase all measurements.Use the \"New image\" button after\n"
                + "loading a new image for resetting all.Using this button will not erase the setted ROIs! \n \n"
                + "---------------------------------------------------------------------------------------\n \n"
                + "written by Joris FA Meys (2009). More info : jorismeys@gmail.com\n"
                + "Canny-deriche and hysteresis plugins written by Thomas Boudier\n"
                + "This plugin is part of the public domain");
    } // End ShowAbout FRAP_Analysis

    private void showPic() {
        ImageJ ij = IJ.getInstance();

        URL url = this.getClass().getResource("/aboutFA.jpg");
        if (url != null) {
            Image img = null;
            try {
                img = ij.createImage((ImageProducer) url.getContent());
            } catch (Exception e) {
            }
            if (img != null) {
                ImagePlus imp = new ImagePlus("", img);
                ImageWindow.centerNextImage();
                imp.show();
            }
        } else {
            showAbout();
        }
    }

    /**
     * Double thresholding
     *
     * @param imageProcessor original image
     * @param highThreshold high threshold
     * @param lowThreshold low threshold
     * @return "trinarised" image
     */
    ImageProcessor trinarise(ByteProcessor imageProcessor, float highThreshold,
            float lowThreshold) {

        int width = imageProcessor.getWidth();
        int height = imageProcessor.getHeight();
        ImageProcessor returnedProcessor = imageProcessor.duplicate();

        for (int x = 0; x < width; x++) {
            for (int y = 0; y < height; y++) {

                float value = returnedProcessor.getPixelValue(x, y);

                if (value >= highThreshold) {
                    returnedProcessor.putPixel(x, y, 255);
                } else if (value >= lowThreshold) {
                    returnedProcessor.putPixel(x, y, 128);
                }
            }
        }
        return returnedProcessor;
    }

    /**
     * Hysteresis thresholding
     *
     * @param imageProcessor original image
     * @return thresholded image
     */
    ImageProcessor hysteresisThresholding(ByteProcessor imageProcessor) {

        int width = imageProcessor.getWidth();
        int height = imageProcessor.getHeight();

        ImageProcessor returnedProcessor = imageProcessor.duplicate();
        boolean change = true;

        while (change) {
            change = false;
            for (int x = 1; x < width - 1; x++) {
                for (int y = 1; y < height - 1; y++) {
                    if (returnedProcessor.getPixelValue(x, y) == MAX_VALUE - 1) {
                        if (returnedProcessor.getPixelValue(x + 1, y) == MAX_VALUE / 2) {
                            change = true;
                            returnedProcessor.putPixelValue(x + 1, y, MAX_VALUE - 1);
                        }
                        if (returnedProcessor.getPixelValue(x - 1, y) == MAX_VALUE / 2) {
                            change = true;
                            returnedProcessor.putPixelValue(x - 1, y, MAX_VALUE - 1);
                        }
                        if (returnedProcessor.getPixelValue(x, y + 1) == MAX_VALUE / 2) {
                            change = true;
                            returnedProcessor.putPixelValue(x, y + 1, MAX_VALUE - 1);
                        }
                        if (returnedProcessor.getPixelValue(x, y - 1) == MAX_VALUE / 2) {
                            change = true;
                            returnedProcessor.putPixelValue(x, y - 1, MAX_VALUE - 1);
                        }
                        if (returnedProcessor.getPixelValue(x + 1, y + 1) == MAX_VALUE / 2) {
                            change = true;
                            returnedProcessor.putPixelValue(x + 1, y + 1, MAX_VALUE - 1);
                        }
                        if (returnedProcessor.getPixelValue(x - 1, y - 1) == MAX_VALUE / 2) {
                            change = true;
                            returnedProcessor.putPixelValue(x - 1, y - 1, MAX_VALUE - 1);
                        }
                        if (returnedProcessor.getPixelValue(x - 1, y + 1) == MAX_VALUE / 2) {
                            change = true;
                            returnedProcessor.putPixelValue(x - 1, y + 1, MAX_VALUE - 1);
                        }
                        if (returnedProcessor.getPixelValue(x + 1, y - 1) == MAX_VALUE / 2) {
                            change = true;
                            returnedProcessor.putPixelValue(x + 1, y - 1, MAX_VALUE - 1);
                        }
                    }
                }
            }
            if (change) {
                for (int x = width - 2; x > 0; x--) {
                    for (int y = height - 2; y > 0; y--) {
                        if (returnedProcessor.getPixelValue(x, y) == MAX_VALUE - 1) {
                            if (returnedProcessor.getPixelValue(x + 1, y) == MAX_VALUE / 2) {
                                change = true;
                                returnedProcessor.putPixelValue(x + 1, y, MAX_VALUE - 1);
                            }
                            if (returnedProcessor.getPixelValue(x - 1, y) == MAX_VALUE / 2) {
                                change = true;
                                returnedProcessor.putPixelValue(x - 1, y, MAX_VALUE - 1);
                            }
                            if (returnedProcessor.getPixelValue(x, y + 1) == MAX_VALUE / 2) {
                                change = true;
                                returnedProcessor.putPixelValue(x, y + 1, MAX_VALUE - 1);
                            }
                            if (returnedProcessor.getPixelValue(x, y - 1) == MAX_VALUE / 2) {
                                change = true;
                                returnedProcessor.putPixelValue(x, y - 1, MAX_VALUE - 1);
                            }
                            if (returnedProcessor.getPixelValue(x + 1, y + 1) == MAX_VALUE / 2) {
                                change = true;
                                returnedProcessor.putPixelValue(x + 1, y + 1, MAX_VALUE - 1);
                            }
                            if (returnedProcessor.getPixelValue(x - 1, y - 1) == MAX_VALUE / 2) {
                                change = true;
                                returnedProcessor.putPixelValue(x - 1, y - 1, MAX_VALUE - 1);
                            }
                            if (returnedProcessor.getPixelValue(x - 1, y + 1) == MAX_VALUE / 2) {
                                change = true;
                                returnedProcessor.putPixelValue(x - 1, y + 1, MAX_VALUE - 1);
                            }
                            if (returnedProcessor.getPixelValue(x + 1, y - 1) == MAX_VALUE / 2) {
                                change = true;
                                returnedProcessor.putPixelValue(x + 1, y - 1, MAX_VALUE - 1);
                            }
                        }
                    }
                }
            }
        }
        // suppression
        for (int x = 0; x < width; x++) {
            for (int y = 0; y < height; y++) {
                if (returnedProcessor.getPixelValue(x, y) == MAX_VALUE / 2) {
                    returnedProcessor.putPixelValue(x, y, MIN_VALUE);
                }
            }
        }
        return returnedProcessor;
    }

    /**
     * Get source tile rectangle.
     *
     * @param x0 X coordinate of the upper left corner point of the target tile
     * rectangle.
     * @param y0 Y coordinate of the upper left corner point of the target tile
     * rectangle.
     * @param w The width of the target tile rectangle.
     * @param h The height of the target tile rectangle.
     * @return The source tile rectangle.
     */
    private Rectangle getSourceTileRectangle(int x0, int y0, int w, int h) {

        int sx0 = x0;
        int sy0 = y0;
        int sw = w;
        int sh = h;

        if (x0 >= halfSizeX) {
            sx0 -= halfSizeX;
            sw += halfSizeX;
        }

        if (y0 >= halfSizeY) {
            sy0 -= halfSizeY;
            sh += halfSizeY;
        }

        if (x0 + w + halfSizeX <= sourceImageWidth) {
            sw += halfSizeX;
        }

        if (y0 + h + halfSizeY <= sourceImageHeight) {
            sh += halfSizeY;
        }

        return new Rectangle(sx0, sy0, sw, sh);
    }

    /**
     * The SPI is used to register this operator in the graph processing
     * framework via the SPI configuration file
     * {@code META-INF/services/org.esa.beam.framework.gpf.OperatorSpi}. This
     * class may also serve as a factory for new operator instances.
     *
     * @see OperatorSpi#createOperator()
     * @see OperatorSpi#createOperator(java.util.Map, java.util.Map)
     */
    public static class Spi extends OperatorSpi {

        public Spi() {
            super(DericheEdgeDetectorOp.class);
            setOperatorUI(DericheEdgeDetectorOpUI.class);
        }
    }
}
