package org.esa.beam.dataio.envisat;

import org.esa.beam.framework.datamodel.MetadataElement;
import org.esa.beam.framework.datamodel.Product;
import org.esa.beam.framework.datamodel.MetadataAttribute;
import org.esa.beam.framework.datamodel.ProductData;

import java.io.File;
import java.util.ArrayList;

/**
  Abstract common metadata from products to be used uniformly by all operators
 */
public final class AsarAbstractMetadata {

    public static final String ABSTRACTED_METADATA_ROOT_NAME = "Abstracted Metadata";

    private final String _productType;
    private final String _version;
    private final File _file;

    AsarAbstractMetadata(String type, String ver, File file) {
        _productType = type;
        _version = ver;
        _file = file;
    }

    /**
     * Abstract common metadata from products to be used uniformly by all operators
     * @param product the product created
     * @param root the product metadata root
     */
    void addAbstractedMetadataHeader(Product product, MetadataElement root) {
        final MetadataElement absRoot = new MetadataElement(ABSTRACTED_METADATA_ROOT_NAME);
        product.getMetadataRoot().addElementAt(absRoot, 0);

        final MetadataElement mph = root.getElement("MPH");
        final MetadataElement sph = root.getElement("SPH");

        final String productType = product.getProductType();
        boolean waveProduct = false;
        if (productType.equals("ASA_WVI_1P") || productType.equals("ASA_WVS_1P") || productType.equals("ASA_WVW_2P")) {
            waveProduct = true;
        }

        MetadataElement mppAds = root.getElement("MAIN_PROCESSING_PARAMS_ADS");
        if(mppAds != null) {
            final MetadataElement ads = mppAds.getElement("MAIN_PROCESSING_PARAMS_ADS.1");
            if(ads != null)
                mppAds = ads;
        }

        // MPH
        addAbstractedAttribute(mph, "PRODUCT", absRoot, "Product name");
        addAbstractedAttribute("PRODUCT_TYPE", _productType, absRoot, "Product type");
        addAbstractedAttribute(sph, "SPH_DESCRIPTOR", absRoot, "Description");
        addAbstractedAttribute("MISSION", getMission(_productType, _file), absRoot, "Satellite mission");
        addAbstractedAttribute("PROC_TIME", mph.getAttributeUTC("PROC_TIME", new ProductData.UTC(0)), absRoot,
                                "Processed time");
        addAbstractedAttribute("Processing system identifier", mph.getAttributeString("SOFTWARE_VER", ""), absRoot,
                                "Processing system identifier");
        addAbstractedAttribute(mph, "CYCLE", absRoot, "Cycle");
        addAbstractedAttribute(mph, "REL_ORBIT", absRoot, "Track");
        addAbstractedAttribute(mph, "ABS_ORBIT", absRoot, "Orbit");
        addAbstractedAttribute("STATE_VECTOR_TIME", mph.getAttributeUTC("STATE_VECTOR_TIME", new ProductData.UTC(0)), absRoot,
                                "Time of orbit state vector");
        addAbstractedAttribute(mph, "VECTOR_SOURCE", absRoot,
                                "State vector source");

        // SPH
        addAbstractedAttribute(sph, "NUM_SLICES", absRoot, "Number of slices");
        if(waveProduct) {
            addAbstractedAttribute("first_line_time", sph.getAttributeUTC("first_cell_time", new ProductData.UTC(0)), absRoot,
                                "First cell time");
            addAbstractedAttribute("last_line_time", sph.getAttributeUTC("last_cell_time", new ProductData.UTC(0)), absRoot,
                                "Last cell time");

            addAbstractedAttribute("SWATH", sph.getAttributeString("SWATH_1", ""), absRoot, "Swath name");
            addAbstractedAttribute(sph, "PASS", absRoot, "ASCENDING or DESCENDING");
            
            String mds1_tx_rx_polar = sph.getAttributeString("tx_rx_polar", "");
            mds1_tx_rx_polar = mds1_tx_rx_polar.replace("/","");
            addAbstractedAttribute("mds1_tx_rx_polar", mds1_tx_rx_polar, absRoot, "Polarization");
            addAbstractedAttribute("mds2_tx_rx_polar", "", absRoot, "Polarization");

            addAbstractedAttribute("srgr_flag", sph.getAttributeInt("SR_GR", 0), "flag", absRoot, "SRGR applied");
            addAbstractedAttribute("ant_elev_corr_flag", sph.getAttributeInt("antenna_corr", 0), "flag", absRoot,
                                "Antenna elevation applied");
            addAbstractedAttribute("map_projection", "", absRoot, "Map projection applied");
        } else {
            addAbstractedAttribute("first_line_time", sph.getAttributeUTC("first_line_time", new ProductData.UTC(0)), absRoot,
                    "First zero doppler azimuth time");
            addAbstractedAttribute("last_line_time", sph.getAttributeUTC("last_line_time", new ProductData.UTC(0)), absRoot,
                    "Last zero doppler azimuth time");

            final double million = 1000000.0;
            addAbstractedAttribute("first_near_lat", sph.getAttributeDouble("first_near_lat", 0) / million, "deg", absRoot, "");
            addAbstractedAttribute("first_near_long", sph.getAttributeDouble("first_near_long", 0) / million, "deg", absRoot, "");
            addAbstractedAttribute("first_far_lat", sph.getAttributeDouble("first_far_lat", 0) / million, "deg", absRoot, "");
            addAbstractedAttribute("first_far_long", sph.getAttributeDouble("first_far_long", 0) / million, "deg", absRoot, "");
            addAbstractedAttribute("last_near_lat", sph.getAttributeDouble("last_near_lat", 0) / million, "deg", absRoot, "");
            addAbstractedAttribute("last_near_long", sph.getAttributeDouble("last_near_long", 0) / million, "deg", absRoot, "");
            addAbstractedAttribute("last_far_lat", sph.getAttributeDouble("last_far_lat", 0) / million, "deg", absRoot, "");
            addAbstractedAttribute("last_far_long", sph.getAttributeDouble("last_far_long", 0) / million, "deg", absRoot, "");

            addAbstractedAttribute(sph, "SWATH", absRoot, "Swath name");
            addAbstractedAttribute(sph, "PASS", absRoot, "ASCENDING or DESCENDING");
            addAbstractedAttribute(sph, "SAMPLE_TYPE", absRoot, "DETECTED or COMPLEX");

            String mds1_tx_rx_polar = sph.getAttributeString("mds1_tx_rx_polar", "");
            mds1_tx_rx_polar = mds1_tx_rx_polar.replace("/","");
            addAbstractedAttribute("mds1_tx_rx_polar", mds1_tx_rx_polar, absRoot, "Polarization");
            String mds2_tx_rx_polar = sph.getAttributeString("mds2_tx_rx_polar", "");
            mds2_tx_rx_polar = mds2_tx_rx_polar.replace("/","");
            addAbstractedAttribute("mds2_tx_rx_polar", mds2_tx_rx_polar, absRoot, "Polarization");
        }

        addAbstractedAttribute(sph, "ALGORITHM", absRoot, "Processing algorithm");
        addAbstractedAttribute(sph, "azimuth_looks", absRoot, "");
        addAbstractedAttribute(sph, "range_looks", absRoot, "");
        addAbstractedAttribute(sph, "range_spacing", absRoot, "Range sample spacing");
        addAbstractedAttribute(sph, "azimuth_spacing", absRoot, "Azimuth sample spacing");

        if(mppAds != null) {
            addAbstractedAttribute("pulse_repetition_frequency", getPulseRepetitionFreq(mppAds), "Hz", absRoot, "PRF");
            addAbstractedAttribute("radar_frequency",
                mppAds.getAttributeDouble("radar_freq", 0) / 1000000.0, "MHz", absRoot, "Radar frequency");
        }
        addAbstractedAttribute(sph, "line_time_interval", absRoot, "");
        addAbstractedAttribute("total_size", (int)(product.getRawStorageSize() / (1024.0f * 1024.0f)), "Mb", absRoot,
                            "Total product size");
        
        //MPP
        if(mppAds != null) {
            addAbstractedAttribute("num_output_lines", product.getSceneRasterHeight(), "", absRoot, "");
            addAbstractedAttribute("num_samples_per_line", product.getSceneRasterWidth(), "", absRoot, "");
            addAbstractedAttribute(mppAds, "srgr_flag", absRoot, "SRGR applied");
            addAbstractedAttribute("avg_scene_height", mppAds.getAttributeDouble("avg_scene_height_ellpsoid", 0),
                    "m", absRoot, "Average ccene height ellipsoid");

            String mapProjection = "";
            String geoRefSystem = "";
            if (productType.contains("APG") || productType.contains("IMG")) {
                final MetadataElement mapGads = root.getElement("MAP_PROJECTION_GADS");
                mapProjection = mapGads.getAttributeString("map_descriptor", "Geocoded");
                geoRefSystem = mapGads.getAttributeString("ellipsoid_name", "WGS84");
            }
            addAbstractedAttribute("map_projection", mapProjection, absRoot, "Map projection applied");

            addAbstractedAttribute("is_terrain_corrected", 0, "flag", absRoot, "orthorectification applied");
            addAbstractedAttribute("dem", "", absRoot, "Digital Elevation Model used");
            addAbstractedAttribute("geo_ref_system", geoRefSystem, absRoot, "geographic reference system");
            addAbstractedAttribute("lat_pixel_res", 0.0, "deg", absRoot, "pixel resolution in geocoded image");
            addAbstractedAttribute("lon_pixel_res", 0.0, "deg", absRoot, "pixel resolution in geocoded image");

            final MetadataElement gg = root.getElement("GEOLOCATION_GRID_ADS");
            if(gg != null) {
                final MetadataElement gg1 = gg.getElement("GEOLOCATION_GRID_ADS.1");
                if(gg1 != null) {
                    final double slantRangeTime = gg1.getAttributeDouble("ASAR_Geo_Grid_ADSR.sd/first_line_tie_points.slant_range_times");
                    final double halfLightSpeed = 299792458.0 / 2.0;
                    final double slantRangeDist = slantRangeTime * halfLightSpeed / 1000000000.0; // slantRangeTime ns to s
                    addAbstractedAttribute("slant_range_to_first_pixel", slantRangeDist, "m", absRoot, "Slant range to 1st data sample");
                }
            }

            addAbstractedAttribute(mppAds, "ant_elev_corr_flag", absRoot, "Antenna elevation applied");
            addAbstractedAttribute(mppAds, "range_spread_comp_flag", absRoot, "range spread compensation applied");
            addAbstractedAttribute("replica_power_corr_flag", ProductData.TYPE_UINT8, "flag",
                    "Replica pulse power correction applied", absRoot);
            addAbstractedAttribute("abs_calibration_flag", ProductData.TYPE_UINT8, "flag",
                    "Product calibrated", absRoot);
            addAbstractedAttribute("calibration_factor",
                mppAds.getAttributeDouble("ASAR_Main_ADSR.sd/calibration_factors.1.ext_cal_fact", 0), "", absRoot,
                    "Calibration constant");
            addAbstractedAttribute("range_sampling_rate",
                mppAds.getAttributeDouble("range_samp_rate", 0) / 1000000.0, "MHz", absRoot, "Range Sampling Rate");
            addAbstractedAttribute("multilook_flag", ProductData.TYPE_UINT8, "flag",
                    "Product multilooked", absRoot);

            addOrbitStateVectors(root, absRoot);
        }

        // add SRGR coefficients if found
        final MetadataElement srgrADS = root.getElement("SR_GR_ADS");
        if(srgrADS != null) {
            addSRGRCoefficients(srgrADS, absRoot);
        }

        final MetadataElement dsd = root.getElement("DSD");
        if(dsd != null) {
            final MetadataElement dsd17 = dsd.getElement("DSD.17");
            if(dsd17 != null) {
                addAbstractedAttribute("external_calibration_file",
                        dsd17.getAttributeString("FILE_NAME", ""), absRoot, "External calibration file used");
            }
            final MetadataElement dsd18 = dsd.getElement("DSD.18");
            if(dsd18 != null) {
                addAbstractedAttribute("orbit_state_vector_file",
                        dsd18.getAttributeString("FILE_NAME", ""), absRoot, "Orbit file used");
            }
        }

    }
    
    public static String getMission(final String productType, final File file) {
        if(productType.startsWith("SAR")) {
            if(file.toString().endsWith("E2"))
                return "ERS2";
            else
                return "ERS1";
        }
        return "ENVISAT";
    }

    /**
     * Adds an attribute from src to dest
     * @param tag the name of the attribute
     * @param value the string value
     * @param dest the destination element
     * @param desc the description
     */
    private static void addAbstractedAttribute(String tag, String value, MetadataElement dest, String desc) {
        if(value == null || value.isEmpty()) 
            value = " ";
        final MetadataAttribute attribute = new MetadataAttribute(tag, ProductData.TYPE_ASCII, 1);
        attribute.getData().setElems(value);
        attribute.setDescription(desc);
        dest.addAttributeFast(attribute);
    }

    /**
     * Adds an attribute from src to dest
     * @param tag the name of the attribute
     * @param value the UTC value
     * @param dest the destination element
     * @param desc the description
     */
    private static void addAbstractedAttribute(String tag, ProductData.UTC value, MetadataElement dest, String desc) {
        final MetadataAttribute attribute = new MetadataAttribute(tag, ProductData.TYPE_UTC, 1);
        attribute.getData().setElems(value.getArray());
        attribute.setUnit("utc");
        attribute.setDescription(desc);
        dest.addAttributeFast(attribute);
    }

    /**
     * Adds an attribute from src to dest
     * @param tag the name of the attribute
     * @param value the UTC value
     * @param unit the unit string
     * @param dest the destination element
     * @param desc the description
     */
    private static void addAbstractedAttribute(String tag, int value, String unit, MetadataElement dest, String desc) {
        final MetadataAttribute attribute = new MetadataAttribute(tag, ProductData.TYPE_INT32, 1);
        attribute.getData().setElemInt(value);
        attribute.setUnit(unit);
        attribute.setDescription(desc);
        dest.addAttributeFast(attribute);
    }

    /**
     * Adds an attribute from src to dest
     * @param tag the name of the attribute
     * @param value the double value
     * @param unit the unit string
     * @param dest the destination element
     * @param desc the description
     */
    private static void addAbstractedAttribute(String tag, double value, String unit, MetadataElement dest, String desc) {
        final MetadataAttribute attribute = new MetadataAttribute(tag, ProductData.TYPE_FLOAT64, 1);
        attribute.getData().setElems( new double[] {value} );
        attribute.setUnit(unit);
        attribute.setDescription(desc);
        dest.addAttributeFast(attribute);
    }

    /**
     * Adds an attribute into dest
     * @param dest the destination element
     * @param tag the name of the attribute
     * @param dataType the ProductData type
     * @param unit The unit
     * @param desc The description
     * @return MetadataAttribute
     */
    private static MetadataAttribute addAbstractedAttribute(String tag, int dataType,
                                               String unit, String desc, MetadataElement dest) {
        final MetadataAttribute attribute = new MetadataAttribute(tag, dataType, 1);
        if(dataType == ProductData.TYPE_ASCII)
            attribute.getData().setElems(" ");
        attribute.setUnit(unit);
        attribute.setDescription(desc);
        attribute.setReadOnly(false);
        dest.addAttributeFast(attribute);
        return attribute;
    }

    /**
     * Adds an attribute from src to dest
     * @param src the source element
     * @param tag the name of the attribute
     * @param dest the destination element
     * @param desc the description
     */
    private static void addAbstractedAttribute(MetadataElement src, String tag, MetadataElement dest, String desc) {
        final MetadataAttribute attrib = src.getAttribute(tag);
        if(attrib != null) {
            MetadataAttribute copiedAttrib = attrib.createDeepClone();
            copiedAttrib.setReadOnly(false);
            copiedAttrib.setDescription(desc);
            dest.addAttributeFast(copiedAttrib);
        }
    }

    private double getPulseRepetitionFreq(MetadataElement mppAds) {
        return mppAds.getAttributeDouble("ASAR_Main_ADSR.sd/image_parameters"+_version+".prf_value", 0);
    }

    private static void addOrbitStateVectors(MetadataElement root, MetadataElement dest) {
        final MetadataElement orbitListElem = new MetadataElement("Orbit_State_Vectors");
        dest.addElement(orbitListElem);

        // get all vectors from all MPP
        final MetadataElement rootMPPAds = root.getElement("MAIN_PROCESSING_PARAMS_ADS");
        final ArrayList<MetadataElement> mppList = new ArrayList<MetadataElement>(10);
        if(rootMPPAds != null) {
            final int numElem = rootMPPAds.getNumElements();
            if(numElem == 0) {
                mppList.add(rootMPPAds);
            } else {
                for(int i=1; i <= numElem; ++i) {
                    final MetadataElement ads = rootMPPAds.getElement("MAIN_PROCESSING_PARAMS_ADS."+i);
                    if(ads != null)
                        mppList.add(ads);
                }
            }
        }

        int i=1;
        for(MetadataElement mppAds : mppList) {
            addVector(mppAds, orbitListElem, "ASAR_Main_ADSR.sd/orbit_state_vectors.1", "orbit_vector"+ i++);
            addVector(mppAds, orbitListElem, "ASAR_Main_ADSR.sd/orbit_state_vectors.2", "orbit_vector"+ i++);
            addVector(mppAds, orbitListElem, "ASAR_Main_ADSR.sd/orbit_state_vectors.3", "orbit_vector"+ i++);
            addVector(mppAds, orbitListElem, "ASAR_Main_ADSR.sd/orbit_state_vectors.4", "orbit_vector"+ i++);
            addVector(mppAds, orbitListElem, "ASAR_Main_ADSR.sd/orbit_state_vectors.5", "orbit_vector"+ i++);
        }
    }

    private static void addVector(MetadataElement mppAds, MetadataElement orbitListElem, String orbitPrefix, String tag) {
        ProductData.UTC utcTime;
        double xPos, yPos, zPos, xVel, yVel, zVel;
        try {
            utcTime = mppAds.getAttributeUTC(orbitPrefix+".state_vect_time_1");
            xPos = mppAds.getAttributeDouble(orbitPrefix+".x_pos_1") / 100.0; // 10^-2 m to m
            yPos = mppAds.getAttributeDouble(orbitPrefix+".y_pos_1") / 100.0; // 10^-2 m to m
            zPos = mppAds.getAttributeDouble(orbitPrefix+".z_pos_1") / 100.0; // 10^-2 m to m
            xVel = mppAds.getAttributeDouble(orbitPrefix+".x_vel_1") / 100000.0; // 10^-5 m/s to m/s
            yVel = mppAds.getAttributeDouble(orbitPrefix+".y_vel_1") / 100000.0; // 10^-5 m/s to m/s
            zVel = mppAds.getAttributeDouble(orbitPrefix+".z_vel_1") / 100000.0; // 10^-5 m/s to m/s
        } catch(Exception e) {
            System.out.println(e.getMessage());
            return;
        }

        final MetadataElement orbitElem = new MetadataElement(tag);
        orbitListElem.addElement(orbitElem);
        addAbstractedAttribute("time", utcTime, orbitElem, "");
        addAbstractedAttribute("x_pos", xPos, "", orbitElem, "");
        addAbstractedAttribute("y_pos", yPos, "", orbitElem, "");
        addAbstractedAttribute("z_pos", zPos, "", orbitElem, "");
        addAbstractedAttribute("x_vel", xVel, "", orbitElem, "");
        addAbstractedAttribute("y_vel", yVel, "", orbitElem, "");
        addAbstractedAttribute("z_vel", zVel, "", orbitElem, "");
    }

    private static void addSRGRCoefficients(MetadataElement srgrAds, MetadataElement dest) {
        final MetadataElement srgrListElem = new MetadataElement("SRGR_Coefficients");
        dest.addElement(srgrListElem);

        int listCnt = 1;
        if(srgrAds.getNumElements() == 0) {
            addCoefficients(srgrAds, srgrListElem, listCnt);
        } else {
            for(MetadataElement srgrSrc : srgrAds.getElements()) {

                addCoefficients(srgrSrc, srgrListElem, listCnt++);
            }
        }
    }

    private static void addCoefficients(MetadataElement srgrSrc, MetadataElement srgrListElem, int listCnt) {
        final ArrayList<Double> coefList = new ArrayList<Double>(5);
        ProductData.UTC utcTime;
        double origin;
        try {
            utcTime = srgrSrc.getAttributeUTC("zero_doppler_time");
            origin = srgrSrc.getAttributeDouble("ground_range_origin");

            final MetadataAttribute srgrCoefAttrib = srgrSrc.getAttribute("srgr_coeff");
            final ProductData data = srgrCoefAttrib.getData();
            for(int i=0; i < data.getNumElems(); ++i) {
                coefList.add(data.getElemDoubleAt(i));
            }
        } catch(Exception e) {
            System.out.println(e.getMessage());
            return;
        }

        final MetadataElement srgrElem = new MetadataElement("srgr_coef_list."+listCnt);
        srgrListElem.addElement(srgrElem);
        addAbstractedAttribute("zero_doppler_time", utcTime, srgrElem, "");
        addAbstractedAttribute("ground_range_origin", origin, "m", srgrElem, "");

        int coefCnt = 1;
        for(Double value : coefList) {
            final MetadataElement coefElem = new MetadataElement("coefficient."+coefCnt);
            ++coefCnt;
            srgrElem.addElement(coefElem);
            addAbstractedAttribute("srgr_coef", value, "", coefElem, "");
        }
    }
}
