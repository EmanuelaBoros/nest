package org.esa.nest.image.processing.segmentation.configuration;

/**
 * Configuration parameters for the ActiveContour plug-in
 *
 * @author Thomas
 * @since 11 May 2004
 * @author Emanuela Boros @updated September 2012
 */
public class ActiveContourConfiguration {

    private double gradThreshold;
    private double maxDisplacement;
    private double maxSearch;
    private double regMin;
    private double regMax;
    private double alphaDeriche;

    /**
     * Constructor for the snakeConfig object
     *
     * @param gt Description of the Parameter
     * @param md Description of the Parameter
     * @param rmin Description of the Parameter
     * @param rmax Description of the Parameter
     * @param ms Description of the Parameter
     * @param alpha Description of the Parameter
     */
    public ActiveContourConfiguration(double gradThreshold, double maxDisplacement,
            double maxSearch, double regMin, double regMax, double alphaDeriche) {
        this.gradThreshold = gradThreshold;
        this.maxDisplacement = maxDisplacement;
        this.maxSearch = maxSearch;
        this.regMin = regMin;
        this.regMax = regMax;
        this.alphaDeriche = alphaDeriche;
    }

    /**
     * Constructor for the ActiveContourConfiguration object
     *
     * @param configuration Description of the Parameter
     */
    public ActiveContourConfiguration(ActiveContourConfiguration configuration) {
        gradThreshold = configuration.getGradThreshold();
        maxDisplacement = configuration.getMaxDisplacement();
        maxSearch = configuration.getMaxSearch();
        regMin = configuration.getRegMin();
        regMax = configuration.getRegMax();
        alphaDeriche = configuration.getAlpha();
    }

    /**
     * Gets the gradThreshold attribute of the ActiveContourConfiguration object
     *
     * @return The gradThreshold value
     */
    public double getGradThreshold() {
        return gradThreshold;
    }

    /**
     * Gets the maxDisplacement attribute of the ActiveContourConfiguration
     * object
     *
     * @return The maxDisplacement value
     */
    public double getMaxDisplacement() {
        return maxDisplacement;
    }

    /**
     * Gets the maxSearch attribute of the ActiveContourConfiguration object
     *
     * @return The maxSearch value
     */
    public double getMaxSearch() {
        return maxSearch;
    }

    /**
     * Gets the regMin attribute of the ActiveContourConfiguration object
     *
     * @return The regMin value
     */
    public double getRegMin() {
        return regMin;
    }

    /**
     * Gets the regMax attribute of the ActiveContourConfiguration object
     *
     * @return The regMax value
     */
    public double getRegMax() {
        return regMax;
    }

    /**
     * Gets the alpha attribute of the ActiveContourConfiguration object
     *
     * @return The alpha value
     */
    public double getAlpha() {
        return alphaDeriche;
    }

    /**
     * Description of the Method
     *
     * @param multiply Description of the Parameter
     */
    public void update(double multiply) {
        alphaDeriche /= multiply;
        //maxDisplacement *= mul;
        //maxSearch *= mul;
        regMax *= multiply;
        regMin *= multiply;
    }
}
