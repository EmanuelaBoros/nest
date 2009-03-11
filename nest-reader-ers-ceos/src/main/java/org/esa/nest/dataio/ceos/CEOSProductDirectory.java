package org.esa.nest.dataio.ceos;

import org.esa.beam.framework.datamodel.*;
import org.esa.beam.framework.dataop.maptransf.Datum;
import org.esa.beam.util.Guardian;
import org.esa.nest.dataio.IllegalBinaryFormatException;
import org.esa.nest.dataio.ReaderUtils;
import org.esa.nest.dataio.ceos.records.BaseRecord;
import org.esa.nest.datamodel.AbstractMetadata;
import org.esa.nest.datamodel.Unit;

import javax.imageio.stream.FileImageInputStream;
import javax.imageio.stream.ImageInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.*;
import java.text.ParseException;

/**

 */
public abstract class CEOSProductDirectory {

    protected boolean isProductSLC = false;
    protected String productType = null;

    protected abstract void readProductDirectory() throws IOException, IllegalBinaryFormatException;

    public abstract Product createProduct() throws IOException, IllegalBinaryFormatException;

    public abstract CEOSImageFile getImageFile(final Band band) throws IOException, IllegalBinaryFormatException;

    public abstract void close() throws IOException;

    public boolean isSLC() {
        return isProductSLC;
    }
    
    public String getSampleType() {
        if(isProductSLC)
            return "COMPLEX";
        else
            return "DETECTED";
    }

    public String getProductType() {
        return productType;
    }

    protected static int getTotalSize(Product product) {
        return (int)(product.getRawStorageSize() / (1024.0f * 1024.0f));
    }

    protected static String getPolarization(String theID) {
        final String id = theID.toUpperCase();
        if(id.contains("HH") || id.contains("H/H") || id.contains("H-H"))
            return "HH";
        else if(id.contains("VV") || id.contains("V/V") || id.contains("V-V"))
            return "VV";
        else if(id.contains("HV") || id.contains("H/V") || id.contains("H-V"))
            return "HV";
        else if(id.contains("VH") || id.contains("V/H") || id.contains("V-H"))
            return "VH";
        return id;
    }

    protected static void addGeoCoding(final Product product, final float[] latCorners, final float[] lonCorners) {

        if(latCorners == null || lonCorners == null) return;

        int gridWidth = 10;
        int gridHeight = 10;

        final float[] fineLatTiePoints = new float[gridWidth*gridHeight];
        ReaderUtils.createFineTiePointGrid(2, 2, gridWidth, gridHeight, latCorners, fineLatTiePoints);

        float subSamplingX = (float)product.getSceneRasterWidth() / (gridWidth - 1);
        float subSamplingY = (float)product.getSceneRasterHeight() / (gridHeight - 1);

        final TiePointGrid latGrid = new TiePointGrid("latitude", gridWidth, gridHeight, 0.5f, 0.5f,
                subSamplingX, subSamplingY, fineLatTiePoints);
        latGrid.setUnit(Unit.DEGREES);

        final float[] fineLonTiePoints = new float[gridWidth*gridHeight];
        ReaderUtils.createFineTiePointGrid(2, 2, gridWidth, gridHeight, lonCorners, fineLonTiePoints);

        final TiePointGrid lonGrid = new TiePointGrid("longitude", gridWidth, gridHeight, 0.5f, 0.5f,
                subSamplingX, subSamplingY, fineLonTiePoints, TiePointGrid.DISCONT_AT_180);
        lonGrid.setUnit(Unit.DEGREES);

        final TiePointGeoCoding tpGeoCoding = new TiePointGeoCoding(latGrid, lonGrid, Datum.WGS_84);

        product.addTiePointGrid(latGrid);
        product.addTiePointGrid(lonGrid);
        product.setGeoCoding(tpGeoCoding);
    }

    protected static ProductData.UTC getProcTime(BaseRecord volDescRec) {
        try {
            final String procDate = volDescRec.getAttributeString("Logical volume preparation date").trim();
            final String procTime = volDescRec.getAttributeString("Logical volume preparation time").trim();

            return ProductData.UTC.parse(procDate + procTime, "yyyyMMddHHmmss");
        } catch(ParseException e) {
            System.out.println(e.toString());
            return new ProductData.UTC(0);
        }
    }

    protected static String getPass(BaseRecord mapProjRec) {
        if(mapProjRec == null) return " ";
        final double heading = mapProjRec.getAttributeDouble("Platform heading at nadir corresponding to scene centre");
        if(heading > 90 && heading < 270) return "DESCENDING";
        else return "ASCENDING";
    }

    protected static ProductData.UTC getUTCScanStartTime(BaseRecord sceneRec) {
        if(sceneRec == null) return new ProductData.UTC(0);
        String startTime = sceneRec.getAttributeString("Zero-doppler azimuth time of first azimuth pixel");
        if(startTime == null || startTime.trim().isEmpty()) {
            startTime = sceneRec.getAttributeString("Scene centre time");
            return AbstractMetadata.parseUTC(startTime.trim(), "yyyymmddHHmmssSSS");
        }
        return AbstractMetadata.parseUTC(startTime);
    }

    protected static ProductData.UTC getUTCScanStopTime(BaseRecord sceneRec) {
        if(sceneRec == null) return new ProductData.UTC(0);
        String endTime = sceneRec.getAttributeString("Zero-doppler azimuth time of last azimuth pixel");
        if(endTime == null || endTime.trim().isEmpty()) {
            endTime = sceneRec.getAttributeString("Scene centre time");
            return AbstractMetadata.parseUTC(endTime.trim(), "yyyymmddHHmmssSSS");
        }
        return AbstractMetadata.parseUTC(endTime);
    }

    protected static double getLineTimeInterval(BaseRecord sceneRec, int sceneHeight) {
        final double startTime = getUTCScanStartTime(sceneRec).getMJD() * 24 * 3600;
        final double stopTime = getUTCScanStopTime(sceneRec).getMJD() * 24 * 3600;
        return (stopTime-startTime) / (sceneHeight-1);
    }

    protected static void addSummaryMetadata(final File summaryFile, final MetadataElement parent) throws IOException {
        if (!summaryFile.exists())
            return;

        final MetadataElement summaryMetadata = new MetadataElement("Summary Information");
        final Properties properties = new Properties();

        properties.load(new FileInputStream(summaryFile));
        final Set unsortedEntries = properties.entrySet();
        final TreeSet sortedEntries = new TreeSet(new Comparator() {
            public int compare(final Object a, final Object b) {
                final Map.Entry entryA = (Map.Entry) a;
                final Map.Entry entryB = (Map.Entry) b;
                return ((String) entryA.getKey()).compareTo((String) entryB.getKey());
            }
        });
        sortedEntries.addAll(unsortedEntries);
        for (Object sortedEntry : sortedEntries) {
            final Map.Entry entry = (Map.Entry) sortedEntry;
            final String data = (String) entry.getValue();
            // strip of double quotes
            final String strippedData = data.substring(1, data.length() - 1);
            final MetadataAttribute attribute = new MetadataAttribute((String) entry.getKey(),
                    new ProductData.ASCII(strippedData), true);
            summaryMetadata.addAttribute(attribute);
        }

        parent.addElement(summaryMetadata);
    }

    protected static void assertSameWidthAndHeightForAllImages(final CEOSImageFile[] imageFiles,
                                                      final int width, final int height) {
        for (int i = 0; i < imageFiles.length; i++) {
            final CEOSImageFile imageFile = imageFiles[i];
            Guardian.assertTrue("_sceneWidth == imageFile[" + i + "].getRasterWidth()",
                                width == imageFile.getRasterWidth());
            Guardian.assertTrue("_sceneHeight == imageFile[" + i + "].getRasterHeight()",
                                height == imageFile.getRasterHeight());
        }
    }

    protected static ImageInputStream createInputStream(final File file) throws IOException {
        return new FileImageInputStream(file);
    }
}