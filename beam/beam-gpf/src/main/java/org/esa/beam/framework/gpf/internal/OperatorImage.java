package org.esa.beam.framework.gpf.internal;

import com.bc.ceres.core.ProgressMonitor;
import org.esa.beam.framework.datamodel.Band;
import org.esa.beam.framework.gpf.Tile;
import org.esa.beam.jai.ImageManager;

import javax.media.jai.ImageLayout;
import javax.media.jai.JAI;
import javax.media.jai.PlanarImage;
import javax.media.jai.SourcelessOpImage;
import java.awt.Rectangle;
import java.awt.image.RenderedImage;
import java.awt.image.WritableRaster;

class OperatorImage extends SourcelessOpImage {

    private static final String DISABLE_TILE_CACHING_PROPERTY = "beam.gpf.disableOperatorTileCaching";

    private final OperatorContext operatorContext;
    private Band targetBand;
    private ProgressMonitor progressMonitor;

    OperatorImage(Band targetBand, OperatorContext operatorContext) {
        this(targetBand, operatorContext, ImageManager.createSingleBandedImageLayout(targetBand));
    }

    private OperatorImage(Band targetBand, OperatorContext operatorContext, ImageLayout imageLayout) {
        super(imageLayout,
              operatorContext.getRenderingHints(),
              imageLayout.getSampleModel(null),
              imageLayout.getMinX(null),
              imageLayout.getMinY(null),
              imageLayout.getWidth(null),
              imageLayout.getHeight(null));
        this.targetBand = targetBand;
        this.operatorContext = operatorContext;
        final boolean disableTileCaching = Boolean.getBoolean(DISABLE_TILE_CACHING_PROPERTY);
        if (disableTileCaching) {
            setTileCache(null);
        } else if (getTileCache() == null) {
            setTileCache(JAI.getDefaultInstance().getTileCache());
        }
    }

    protected OperatorContext getOperatorContext() {
        return operatorContext;
    }

    protected Band getTargetBand() {
        return targetBand;
    }


    @Override
    protected void computeRect(PlanarImage[] ignored, WritableRaster tile, Rectangle destRect) {

        long nanos1 = System.nanoTime();

        Tile targetTile;
        if (getOperatorContext().isComputing(getTargetBand())) {
            targetTile = createTargetTile(getTargetBand(), tile, destRect);
        } else if (requiresAllBands()) {
            targetTile = getOperatorContext().getSourceTile(getTargetBand(), destRect, getProgressMonitor());
        } else {
            targetTile = null;
        }
        // computeTile() may have been deactivated
        //if (operatorContext.isComputeTileMethodUsable()) {
		if (targetTile != null && getOperatorContext().isComputeTileMethodUsable()) {
            getOperatorContext().getOperator().computeTile(getTargetBand(), targetTile, getProgressMonitor());
        }

        long nanos2 = System.nanoTime();
        updatePerformanceMetrics(nanos1, nanos2, destRect);
    }

    protected boolean requiresAllBands() {
        return operatorContext.requiresAllBands();
    }

    protected void updatePerformanceMetrics(long nanos1, long nanos2, Rectangle destRect) {
        double targetNanosPerPixel = (double) (nanos2 - nanos1) / (double) (destRect.width * destRect.height);
        operatorContext.getPerformanceMetric().updateTarget(targetNanosPerPixel);

        double sourceNanosPerPixel = operatorContext.getSourceNanosPerPixel();
        operatorContext.getPerformanceMetric().updateSource(sourceNanosPerPixel);
    }

    protected ProgressMonitor getProgressMonitor() {
        return progressMonitor != null ? progressMonitor : ProgressMonitor.NULL;
    }

    protected void setProgressMonitor(ProgressMonitor pm) {
        this.progressMonitor = pm;
    }

    @Override
    public synchronized void dispose() {
        //targetBand = null;
        progressMonitor = null;
        super.dispose();
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder(getClass().getSimpleName());
        sb.append("[");
        sb.append(operatorContext.getOperatorSpi().getOperatorAlias());
        sb.append(",");
        if (targetBand != null) {
            sb.append(targetBand.getName());
        }
        sb.append("]");
        return sb.toString();
    }

    static ProgressMonitor setProgressMonitor(RenderedImage image, ProgressMonitor pm) {
        ProgressMonitor oldPm = ProgressMonitor.NULL;
        if (image instanceof OperatorImage) {
            final OperatorImage opImage = (OperatorImage) image;
            oldPm = opImage.getProgressMonitor();
            opImage.setProgressMonitor(pm);
        }
        return oldPm;
    }

    protected static TileImpl createTargetTile(Band band, WritableRaster targetTileRaster, Rectangle targetRectangle) {
        return new TileImpl(band, targetTileRaster, targetRectangle, true);
    }
}
