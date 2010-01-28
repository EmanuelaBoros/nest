package org.esa.beam.visat.toolviews.layermanager.layersrc.shapefile;

import com.vividsolutions.jts.geom.Geometry;

import org.esa.beam.framework.datamodel.Product;
import org.esa.beam.framework.datamodel.RasterDataNode;
import org.esa.beam.jai.ImageManager;
import org.esa.beam.util.FeatureCollectionClipper;
import org.geotools.data.DataStore;
import org.geotools.data.DataStoreFinder;
import org.geotools.data.FeatureSource;
import org.geotools.data.shapefile.ShapefileDataStoreFactory;
import org.geotools.factory.CommonFactoryFinder;
import org.geotools.feature.FeatureCollection;
import org.geotools.styling.SLDParser;
import org.geotools.styling.Style;
import org.opengis.feature.simple.SimpleFeature;
import org.opengis.feature.simple.SimpleFeatureType;
import org.opengis.referencing.crs.CoordinateReferenceSystem;

import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.util.HashMap;
import java.util.Map;

/**
 * Unstable API. Use at own risk.
 */
public class ShapefileUtils {
    
    private static final org.geotools.styling.StyleFactory styleFactory = CommonFactoryFinder.getStyleFactory(null);

    /**
     * Loads a shapefile into a feature collection. The shapefile is clipped to the geometry of the given raster.
     *
     * @param file         The shapefile.
     * @param targetRaster A geocoded raster.
     *
     * @return The shapefile as a feature collection clipped to the geometry of the raster.
     *
     * @throws IOException if the shapefile could not be read.
     */
    public static FeatureCollection<SimpleFeatureType, SimpleFeature> loadShapefile(File file,
                                                                                    RasterDataNode targetRaster) throws IOException {
        return loadShapefile(file.toURI().toURL(), targetRaster);
    }

    /**
     * Loads a shapefile into a feature collection. The shapefile is clipped to the geometry of the given product.
     *
     * @param file         The shapefile.
     * @param product      A geocoded product.
     *
     * @return The shapefile as a feature collection clipped to the geometry of the product.
     *
     * @throws IOException if the shapefile could not be read.
     */
    public static FeatureCollection<SimpleFeatureType, SimpleFeature> loadShapefile(File file,
                                                                                    Product product) throws IOException {
        final URL url = file.toURI().toURL();
        final CoordinateReferenceSystem targetCrs = ImageManager.getModelCrs(product.getGeoCoding());
        final Geometry clipGeometry = FeatureCollectionClipper.createGeoBoundaryPolygon(product);
        return createFeatureCollection(url, targetCrs, clipGeometry);
    }
    
    /**
     * Loads a shapefile into a feature collection. The shapefile is clipped to the geometry of the given raster.
     *
     * @param url      The URL of the shapefile.
     * @param targetRaster A geocoded raster.
     *
     * @return The shapefile as a feature collection clipped to the geometry of the raster.
     *
     * @throws IOException if the shapefile could not be read.
     */
    public static FeatureCollection<SimpleFeatureType, SimpleFeature> loadShapefile(URL url,
                                                                                    RasterDataNode targetRaster) throws IOException {
        final CoordinateReferenceSystem targetCrs = ImageManager.getModelCrs(targetRaster.getGeoCoding());
        final Geometry clipGeometry = FeatureCollectionClipper.createGeoBoundaryPolygon(targetRaster.getProduct());
        return createFeatureCollection(url, targetCrs, clipGeometry);
    }

    public static FeatureCollection<SimpleFeatureType, SimpleFeature> createFeatureCollection(URL url,
                                                                                              CoordinateReferenceSystem targetCrs,
                                                                                              Geometry clipGeometry
    ) throws IOException {
        FeatureCollection<SimpleFeatureType, SimpleFeature> featureCollection;
        FeatureSource<SimpleFeatureType, SimpleFeature> featureSource = getFeatureSource(url);
        featureCollection = featureSource.getFeatures();
        featureCollection = FeatureCollectionClipper.doOperation(featureCollection, clipGeometry, null, targetCrs);
        return featureCollection;
    }

    public static FeatureSource<SimpleFeatureType, SimpleFeature> getFeatureSource(URL url) throws IOException {
        Map<String, Object> map = new HashMap<String, Object>();
        map.put(ShapefileDataStoreFactory.URLP.key, url);
        map.put(ShapefileDataStoreFactory.CREATE_SPATIAL_INDEX.key, Boolean.TRUE);
        DataStore shapefileStore = DataStoreFinder.getDataStore(map);
        String typeName = shapefileStore.getTypeNames()[0]; // Shape files do only have one type name
        FeatureSource<SimpleFeatureType, SimpleFeature> featureSource;
        featureSource = shapefileStore.getFeatureSource(typeName);
        return featureSource;
    }
    
    /**
     * Figure out the URL for the "sld" file
     */
    private static File getSLDFile(File shapeFile) {
        String filename = shapeFile.getAbsolutePath();
        if (filename.endsWith(".shp") || filename.endsWith(".dbf")
            || filename.endsWith(".shx")) {
            filename = filename.substring(0, filename.length() - 4);
            filename += ".sld";
        } else if (filename.endsWith(".SHP") || filename.endsWith(".DBF")
                   || filename.endsWith(".SHX")) {
            filename = filename.substring(0, filename.length() - 4);
            filename += ".SLD";
        }
        return new File(filename);
    }
    
    public static Style[] loadSLD(File shapeFile) {
        File sld = ShapefileUtils.getSLDFile(shapeFile);
        if (sld.exists()) {
            return createFromSLD(sld);
        } else {
            return new Style[0];
        }
    }
    
    private static Style[] createFromSLD(File sld) {
        try {
            SLDParser stylereader = new SLDParser(styleFactory, sld.toURI().toURL());
            return stylereader.readXML();
        } catch (IOException e) {
            e.printStackTrace();
        }
        return new Style[0];
    }
}
