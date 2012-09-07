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
import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.gui.Roi;
import ij.io.FileSaver;
import ij.plugin.frame.RoiManager;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;
import java.awt.Color;
import java.awt.Rectangle;
import java.awt.image.BufferedImage;
import java.awt.image.RenderedImage;
import java.util.HashMap;
import java.util.Map;
import org.esa.beam.framework.datamodel.Band;
import org.esa.beam.framework.datamodel.Product;
import org.esa.beam.framework.datamodel.ProductData;
import org.esa.beam.framework.gpf.Operator;
import org.esa.beam.framework.gpf.OperatorException;
import org.esa.beam.framework.gpf.OperatorSpi;
import org.esa.beam.framework.gpf.Tile;
import org.esa.beam.framework.gpf.annotations.OperatorMetadata;
import org.esa.beam.framework.gpf.annotations.Parameter;
import org.esa.beam.framework.gpf.annotations.SourceProduct;
import org.esa.beam.framework.gpf.annotations.TargetProduct;
import org.esa.nest.gpf.OperatorUtils;
import org.esa.nest.image.processing.utils.SnakeConfig;
import org.esa.nest.image.processing.utils.SnakeConfigDriver;

@OperatorMetadata(alias = "ActiveContour",
category = "SAR Tools\\Image Processing",
description = "ActiveContour")
public class ActiveContourOp extends Operator {

    public static float[] probabilityHistogram;
    final static int MAX_VALUE = 0;
    final static int MIN_VALUE = 256;
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
    private static ByteProcessor fullByteProcessor;
    /**
     *
     */
    Roi originalROI = null;
    Color colorDraw = Color.white;
    /**
     *
     */
    SnakeConfigDriver configDriver;
    @Parameter(description = "Number of Iterations", defaultValue = "100", label = "nIterations")
    int nIterations = 50;
    // step to display snake
    int step = 1;
    // threshold of edges
    @Parameter(description = "Gradient threshold", defaultValue = "100", label = "GradientThreshold")
    private int gradientThreshold = 5;
    // how far to look for edges
    int DistMAX = Prefs.getInt("ABSnake_DistSearch.int", 100);
    // maximum displacement
    double force = 5.0;
    // regularization factors, min and max
    double dRegularization = 5.0;
    double dMinRegularization, dMaxRegularization;
    // first and last slice to process
    int slice1, slice2;
    // misc options
    boolean showgrad = false;
    boolean savecoords = false;
    boolean createsegimage = false;
    boolean advanced = false;
    boolean propagate = true;
    boolean movie = false;

    /**
     * Initializes this operator and sets the one and only target product.
     * <p>The target product can be either defined by a field of type {@link org.esa.beam.framework.datamodel.Product}
     * annotated with the
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

            initializeActiveContoursProcessing(sourceBand, sourceRaster, targetTile, x0, y0, w, h, pm);

        } catch (Throwable e) {
            OperatorUtils.catchOperatorException(getId(), e);
        } finally {
            pm.done();
        }
    }

    /**
     * Apply Otsu Thresholding
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
    private synchronized void initializeActiveContoursProcessing(final Band sourceBand, final Tile sourceRaster,
            final Tile targetTile, final int x0, final int y0, final int w, final int h,
            final ProgressMonitor pm) {

        if (!processed) {
            final RenderedImage fullRenderedImage = sourceBand.getSourceImage().getImage(0);
            final BufferedImage fullBufferedImage = new BufferedImage(sourceBand.getSceneRasterWidth(),
                    sourceBand.getSceneRasterHeight(),
                    BufferedImage.TYPE_USHORT_GRAY);
            fullBufferedImage.setData(fullRenderedImage.getData());

            fullImagePlus = new ImagePlus(sourceBand.getDisplayName(), fullBufferedImage);

            final ImageProcessor fullImageProcessor = fullImagePlus.getProcessor();

            fullByteProcessor = (ByteProcessor) fullImageProcessor.convertToByte(true);

            showDialog();

            RoiManager roiManager = RoiManager.getInstance();

            if (roiManager == null) {
                roiManager = new RoiManager();
                roiManager.setVisible(true);
                originalROI = fullImagePlus.getRoi();
                if (originalROI == null) {
                    IJ.showMessage("Roi required");
                } else {
                    roiManager.add(fullImagePlus, originalROI, 0);
                }
            }

            processed = true;
        }
        configDriver = new SnakeConfigDriver();
        setAdvancedParameters();

        final Rectangle srcTileRectangle = sourceRaster.getRectangle();
        Roi currentROI = new Roi(sourceRaster.getRectangle());

        ImageProcessor aPartProcessor = fullByteProcessor.duplicate();

        aPartProcessor.setRoi(srcTileRectangle);

        ImageProcessor roiImageProcessor = aPartProcessor.crop();

        configDriver = new SnakeConfigDriver();
        setAdvancedParameters();

        dMinRegularization = dRegularization / 2.0;
        dMaxRegularization = dRegularization;

        ActiveContour currentActivecontour = processActiveContour(fullImagePlus,
                currentROI, roiImageProcessor, w,
                x0 + "," + y0);

        currentActivecontour.drawContour(aPartProcessor, colorDraw, 2);

        final ProductData trgData = targetTile.getDataBuffer();
        final ProductData sourceData = ProductData.createInstance((byte[]) roiImageProcessor.getPixels());

        final int maxY = y0 + h;
        final int maxX = x0 + w;
        for (int y = y0; y < maxY; ++y) {
            for (int x = x0; x < maxX; ++x) {

                trgData.setElemFloatAt(targetTile.getDataBufferIndex(x, y),
                        sourceData.getElemFloatAt(sourceRaster.getDataBufferIndex(x, y)));
            }
        }
    }

    private boolean showDialog() {
        String[] colors = {"Red", "Green", "Blue", "Cyan", "Magenta", "Yellow", "Black", "White"};
        int indexcol = 0;
        GenericDialog gd = new GenericDialog("Snake");
        gd.addNumericField("Gradient_threshold:", gradientThreshold, 0);
        gd.addNumericField("Number_of_iterations:", nIterations, 0);
        gd.addNumericField("Step_result_show:", step, 0);
        //if (profondeur == 1) {
        gd.addCheckbox("Save intermediate images", movie);
        //}
//        if (profondeur > 1) {
//            gd.addNumericField("First_slice:", slice1, 0);
//            gd.addNumericField("Last_slice:", slice2, 0);
//            gd.addCheckbox("Propagate roi", propagate);
//        }
        gd.addChoice("Draw_color:", colors, colors[indexcol]);
        gd.addCheckbox("Save_coords", savecoords);
        gd.addCheckbox("Create_seg_image", createsegimage);
        //gd.addCheckbox("Advanced_options", advanced);
        // show dialog
        gd.showDialog();

        // threshold of edge
        gradientThreshold = (int) gd.getNextNumber();

        // number of iterations
        nIterations = (int) gd.getNextNumber();
        // step of display
        step = (int) gd.getNextNumber();
        //if (profondeur == 1) {
        movie = gd.getNextBoolean();
        //}
        if (step > nIterations - 1) {
            IJ.showStatus("Warning : show step too big\n\t step assignation 1");
            step = 1;
        }
//        if (profondeur > 1) {
//            slice1 = (int) gd.getNextNumber();
//            slice2 = (int) gd.getNextNumber();
//            propagate = gd.getNextBoolean();
//        }
        // color choice of display
        indexcol = gd.getNextChoiceIndex();
        switch (indexcol) {
            case 0:
                colorDraw = Color.red;
                break;
            case 1:
                colorDraw = Color.green;
                break;
            case 2:
                colorDraw = Color.blue;
                break;
            case 3:
                colorDraw = Color.cyan;
                break;
            case 4:
                colorDraw = Color.magenta;
                break;
            case 5:
                colorDraw = Color.yellow;
                break;
            case 6:
                colorDraw = Color.black;
                break;
            case 7:
                colorDraw = Color.white;
                break;
            default:
                colorDraw = Color.yellow;
        }
        savecoords = gd.getNextBoolean();
        createsegimage = gd.getNextBoolean();
        //advanced = gd.getNextBoolean();

        return !gd.wasCanceled();
    }

    /**
     * Main processing method for the Snake_deriche_ object
     *
     * @param imageProcessor image
     */
//    public void letOffActiveContours() {
//
//        stack = fullImagePlus.getStack();
//
//        stackSize = stack.getSize();
//        stackWidth = stack.getWidth();
//        stackHeight = stack.getHeight();
//
//        slice1 = 1;
//        slice2 = stackSize;
//
//        RoiManager roiManager = RoiManager.getInstance();
//        if (roiManager == null) {
//            roiManager = new RoiManager();
//            roiManager.setVisible(true);
//            originalROI = fullImagePlus.getRoi();
//            if (originalROI == null) {
//                IJ.showMessage("Roi required");
//            } else {
//                roiManager.add(fullImagePlus, originalROI, 0);
//            }
//        }
//
//        nROI = roiManager.getCount();
//
//        final Roi[] originalROIs = roiManager.getRoisAsArray();
//        final Roi[] currentROIs = new Roi[nROI];
//
//        Roi[] RoisResult = new Roi[nROI];
//
//        System.arraycopy(originalROIs, 0, currentROIs, 0, nROI);
//
//        configDriver = new SnakeConfigDriver();
//        setAdvancedParameters();
//
//        dMinRegularization = dRegularization / 2.0;
//        dMaxRegularization = dRegularization;
//        // ?
//        // init result
//        resultStack = new ImageStack(stackWidth, stackHeight, java.awt.image.ColorModel.getRGBdefault());
//        if (createsegimage) {
//            pile_seg = new ImageStack(stackWidth, stackHeight);
//        }
//
//        int nbcpu = 1;
//        Thread[] threads = new Thread[nbcpu];
//        final AtomicInteger k = new AtomicInteger(0);
//        final ActiveContour[] activeContours = new ActiveContour[originalROIs.length];
//
//        // for all slices
//
//        // display in RGB color
//        final ColorProcessor[] colorProcessors = new ColorProcessor[originalROIs.length];
//        final ImagePlus[] imagesPlus = new ImagePlus[originalROIs.length];
//
//        int iDirection = slice1 < slice2 ? 1 : -1;
//        for (int z = slice1; z != (slice2 + iDirection); z += iDirection) {
//            final int zz = z;
//            k.set(0);
//
//            for (int i = 0; i < originalROIs.length; i++) {
//                colorProcessors[i] = (ColorProcessor) (resultStack.getProcessor(zz).duplicate());
//                imagesPlus[i] = new ImagePlus("Roi " + i, colorProcessors[i]);
//            }
//            for (int t = 0; t < threads.length; t++) {
//                threads[t] = new Thread() {
//
//                    @Override
//                    public void run() {
//                        IJ.wait(1000);
//                        Roi roi = null;
//                        for (int i = k.getAndIncrement(); i < originalROIs.length; i = k.getAndIncrement()) {
//
//                            if (propagate) {
//                                roi = currentROIs[i];
//                            } else {
//                                roi = originalROIs[i];
//                            }
//                            IJ.log("processing slice " + zz + " with roi " + i);
//
//                            activeContours[i] = processSnake(imagesPlus[i], roi, zz, i + 1);
//                        }
//                    }
//                };
//            }
//
//            for (int ithread = 0; ithread < threads.length; ++ithread) {
//                threads[ithread].setPriority(Thread.NORM_PRIORITY);
//                threads[ithread].start();
//            }
//            try {
//                for (int ithread = 0; ithread < threads.length; ++ithread) {
//                    threads[ithread].join();
//                }
//            } catch (InterruptedException ie) {
//                throw new RuntimeException(ie);
//            }
//
//            RoiEncoder saveRoi;
//            ColorProcessor imageDraw = (ColorProcessor) (resultStack.getProcessor(zz).duplicate());
//            for (int i = 0; i < originalROIs.length; i++) {
//                activeContours[i].drawContour(imageDraw, colorDraw, 1);
//                imagesPlus[i].hide();
//                RoisResult[i] = activeContours[i].createRoi();
//                RoisResult[i].setName("res-" + i);
//                currentROIs[i] = activeContours[i].createRoi();
//            }
//            resultStack.setPixels(imageDraw.getPixels(), z);
//
//            if (createsegimage) {
//                ByteProcessor seg = new ByteProcessor(pile_seg.getWidth(), pile_seg.getHeight());
//                ByteProcessor tmp;
//                for (int i = 0; i < originalROIs.length; i++) {
//                    tmp = activeContours[i].segmentation(seg.getWidth(), seg.getHeight(), i + 1);
//                    seg.copyBits(tmp, 0, 0, Blitter.ADD);
//                }
//                seg.resetMinAndMax();
//                pile_seg.addSlice("Seg " + z, seg);
//            }
//
////            if (savecoords) {
////                for (int i = 0; i < RoisOrig.length; i++) {
////                    try {
////                        activeContours[i].writeCoordinates("ABSnake-r" + (i + 1) + "-z", zz, resXY);
////                        saveRoi = new RoiEncoder("ABSnake-r" + (i + 1) + "-z" + zz + ".roi");
////                        saveRoi.write(RoisResult[i]);
////                    } catch (IOException ex) {
////                        Logger.getLogger(ActiveContour.class.getName()).log(Level.SEVERE, null, ex);
////                    }
////                }
////            } 
//        }
////        new ImagePlus("Draw", resultStack).show();
////        if (createsegimage) {
////            new ImagePlus("Seg", pile_seg).show();
////        }
//    }
    /**
     * Dialog advanced
     *
     * @return dialog ok ?
     */
    private void setAdvancedParameters() {
        configDriver.setMaxDisplacement(Prefs.get("ABSnake_DisplMin.double", 0.1), Prefs.get("ABSnake_DisplMax.double", 2.0));
        configDriver.setInvAlphaD(Prefs.get("ABSnake_InvAlphaMin.double", 0.5), Prefs.get("ABSnake_InvAlphaMax.double", 2.0));
        configDriver.setReg(Prefs.get("ABSnake_RegMin.double", 0.1), Prefs.get("ABSnake_RegMax.double", 2.0));
        configDriver.setStep(Prefs.get("ABSnake_MulFactor.double", 0.99));
    }

    /**
     * do the snake algorithm on all images
     *
     * @param image RGB image to display the snake
     * @param numSlice which image of the stack
     */
    public ActiveContour processActiveContour(ImagePlus imagePlus, Roi currentROI,
            ImageProcessor roiImageProcessor, int numSlice, String numRoi) {

        int i;

        SnakeConfig config;

        ActiveContour activeContour = new ActiveContour();
        activeContour.Init(currentROI);
        activeContour.setOriginalImage(roiImageProcessor);

        IJ.showStatus("Calculating snake...");
        if (step > 0) {
            imagePlus.show();
        }

        double InvAlphaD = configDriver.getInvAlphaD(false);
        double regMax = configDriver.getReg(false);
        double regMin = configDriver.getReg(true);
        double DisplMax = configDriver.getMaxDisplacement(false);
        double mul = configDriver.getStep();

        config = new SnakeConfig(gradientThreshold, DisplMax, DistMAX, regMin, regMax, 1.0 / InvAlphaD);
        activeContour.setConfig(config);
        // compute image gradient
        activeContour.computeGrad(roiImageProcessor);

        IJ.resetEscape();

        FileSaver fs = new FileSaver(imagePlus);

        double dist0 = 0.0;
        double dist;
        //double InvAlphaD0 = InvAlphaD;
        for (i = 0; i < nIterations; i++) {
            if (IJ.escapePressed()) {
                break;
            }
            dist = activeContour.process();
            if ((dist >= dist0) && (dist < force)) {
                activeContour.computeGrad(roiImageProcessor);
                config.update(mul);
            }
            dist0 = dist;

            if ((step > 0) && ((i % step) == 0)) {
                IJ.showStatus("Show intermediate result (iteration n" + (i + 1) + ")");
                ByteProcessor image2 = (ByteProcessor) (roiImageProcessor.duplicate());

                activeContour.drawContour(image2, colorDraw, 1);
                imagePlus.setProcessor("", image2);
                imagePlus.setTitle(fullImagePlus.getTitle() + " roi " + numRoi
                        + " (iteration n" + (i + 1) + ")");
                imagePlus.updateAndRepaintWindow();
                if (movie) {
                    fs = new FileSaver(imagePlus);
                    fs.saveAsTiff("ABsnake-t" + i + "-r" + numRoi + "-z"
                            + numSlice + ".tif");
                }
            }
        }
        return activeContour;
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
                    if (returnedProcessor.getPixelValue(x, y) == 255) {
                        if (returnedProcessor.getPixelValue(x + 1, y) == 128) {
                            change = true;
                            returnedProcessor.putPixelValue(x + 1, y, 255);
                        }
                        if (returnedProcessor.getPixelValue(x - 1, y) == 128) {
                            change = true;
                            returnedProcessor.putPixelValue(x - 1, y, 255);
                        }
                        if (returnedProcessor.getPixelValue(x, y + 1) == 128) {
                            change = true;
                            returnedProcessor.putPixelValue(x, y + 1, 255);
                        }
                        if (returnedProcessor.getPixelValue(x, y - 1) == 128) {
                            change = true;
                            returnedProcessor.putPixelValue(x, y - 1, 255);
                        }
                        if (returnedProcessor.getPixelValue(x + 1, y + 1) == 128) {
                            change = true;
                            returnedProcessor.putPixelValue(x + 1, y + 1, 255);
                        }
                        if (returnedProcessor.getPixelValue(x - 1, y - 1) == 128) {
                            change = true;
                            returnedProcessor.putPixelValue(x - 1, y - 1, 255);
                        }
                        if (returnedProcessor.getPixelValue(x - 1, y + 1) == 128) {
                            change = true;
                            returnedProcessor.putPixelValue(x - 1, y + 1, 255);
                        }
                        if (returnedProcessor.getPixelValue(x + 1, y - 1) == 128) {
                            change = true;
                            returnedProcessor.putPixelValue(x + 1, y - 1, 255);
                        }
                    }
                }
            }
            if (change) {
                for (int x = width - 2; x > 0; x--) {
                    for (int y = height - 2; y > 0; y--) {
                        if (returnedProcessor.getPixelValue(x, y) == 255) {
                            if (returnedProcessor.getPixelValue(x + 1, y) == 128) {
                                change = true;
                                returnedProcessor.putPixelValue(x + 1, y, 255);
                            }
                            if (returnedProcessor.getPixelValue(x - 1, y) == 128) {
                                change = true;
                                returnedProcessor.putPixelValue(x - 1, y, 255);
                            }
                            if (returnedProcessor.getPixelValue(x, y + 1) == 128) {
                                change = true;
                                returnedProcessor.putPixelValue(x, y + 1, 255);
                            }
                            if (returnedProcessor.getPixelValue(x, y - 1) == 128) {
                                change = true;
                                returnedProcessor.putPixelValue(x, y - 1, 255);
                            }
                            if (returnedProcessor.getPixelValue(x + 1, y + 1) == 128) {
                                change = true;
                                returnedProcessor.putPixelValue(x + 1, y + 1, 255);
                            }
                            if (returnedProcessor.getPixelValue(x - 1, y - 1) == 128) {
                                change = true;
                                returnedProcessor.putPixelValue(x - 1, y - 1, 255);
                            }
                            if (returnedProcessor.getPixelValue(x - 1, y + 1) == 128) {
                                change = true;
                                returnedProcessor.putPixelValue(x - 1, y + 1, 255);
                            }
                            if (returnedProcessor.getPixelValue(x + 1, y - 1) == 128) {
                                change = true;
                                returnedProcessor.putPixelValue(x + 1, y - 1, 255);
                            }
                        }
                    }
                }
            }
        }
        // suppression
        for (int x = 0; x < width; x++) {
            for (int y = 0; y < height; y++) {
                if (returnedProcessor.getPixelValue(x, y) == 128) {
                    returnedProcessor.putPixelValue(x, y, 0);
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
     * @param w The stackWidth of the target tile rectangle.
     * @param h The stackHeight of the target tile rectangle.
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
            super(ActiveContourOp.class);
            setOperatorUI(ActiveContourOpUI.class);
        }
    }
}
