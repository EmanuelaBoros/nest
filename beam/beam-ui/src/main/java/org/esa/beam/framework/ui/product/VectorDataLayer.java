/*
 * $Id: VectorDataLayer.java,v 1.4 2009-12-21 16:13:40 lveci Exp $
 *
 * Copyright (C) 2008 by Brockmann Consult (info@brockmann-consult.de)
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation. This program is distributed in the hope it will
 * be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */
package org.esa.beam.framework.ui.product;

import com.bc.ceres.binding.PropertySet;
import com.bc.ceres.glayer.Layer;
import com.bc.ceres.glayer.LayerContext;
import com.bc.ceres.grender.Rendering;
import com.bc.ceres.swing.figure.AbstractFigureChangeListener;
import com.bc.ceres.swing.figure.Figure;
import com.bc.ceres.swing.figure.FigureChangeEvent;
import com.bc.ceres.swing.figure.FigureCollection;
import com.bc.ceres.swing.figure.support.DefaultFigureCollection;
import org.esa.beam.framework.datamodel.ProductNode;
import org.esa.beam.framework.datamodel.ProductNodeEvent;
import org.esa.beam.framework.datamodel.ProductNodeListenerAdapter;
import org.esa.beam.framework.datamodel.VectorDataNode;
import org.geotools.feature.FeatureCollection;
import org.geotools.feature.FeatureIterator;
import org.opengis.feature.simple.SimpleFeature;
import org.opengis.feature.simple.SimpleFeatureType;

import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.geom.Rectangle2D;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;

public class VectorDataLayer extends Layer {

    private static final VectorDataLayerType TYPE = new VectorDataLayerType();
    private VectorDataNode vectorDataNode;
    private final SimpleFeatureFigureFactory figureFactory;
    private FigureCollection figureCollection;
    private VectorDataChangeHandler vectorDataChangeHandler;
    private boolean adjustingFigures;


    public VectorDataLayer(LayerContext ctx, VectorDataNode vectorDataNode) {
        this(TYPE, vectorDataNode, TYPE.createLayerConfig(ctx));
        getConfiguration().setValue(VectorDataLayerType.PROPERTY_NAME_VECTOR_DATA, vectorDataNode.getName());
    }

    VectorDataLayer(VectorDataLayerType vectorDataLayerType, VectorDataNode vectorDataNode, PropertySet configuration) {
        super(vectorDataLayerType, configuration);

        this.vectorDataNode = vectorDataNode;
        setName(vectorDataNode.getName());
        figureFactory = new SimpleFeatureFigureFactory(vectorDataNode.getFeatureCollection());
        figureCollection = new DefaultFigureCollection();
        updateFigureCollection();

        vectorDataChangeHandler = new VectorDataChangeHandler();
        vectorDataNode.getProduct().addProductNodeListener(vectorDataChangeHandler);
        figureCollection.addChangeListener(new FigureCollectionChangeListener());
    }

    public VectorDataNode getVectorDataNode() {
        return vectorDataNode;
    }

    @Override
    protected void disposeLayer() {
        vectorDataNode.getProduct().removeProductNodeListener(vectorDataChangeHandler);
        vectorDataNode = null;
        super.disposeLayer();
    }

    private void updateFigureCollection() {
        FeatureCollection<SimpleFeatureType, SimpleFeature> featureCollection = vectorDataNode.getFeatureCollection();

        Figure[] figures = figureCollection.getFigures();
        Map<SimpleFeature, SimpleFeatureFigure> figureMap = new HashMap<SimpleFeature, SimpleFeatureFigure>();
        for (Figure figure : figures) {
            if (figure instanceof SimpleFeatureFigure) {
                SimpleFeatureFigure simpleFeatureFigure = (SimpleFeatureFigure) figure;
                figureMap.put(simpleFeatureFigure.getSimpleFeature(), simpleFeatureFigure);
            }
        }

        FeatureIterator<SimpleFeature> featureIterator = featureCollection.features();
        while (featureIterator.hasNext()) {
            SimpleFeature simpleFeature = featureIterator.next();
            SimpleFeatureFigure featureFigure = figureMap.get(simpleFeature);
            if (featureFigure != null) {
                figureMap.remove(simpleFeature);
            } else {
                featureFigure = getFigureFactory().createSimpleFeatureFigure(simpleFeature, vectorDataNode.getDefaultCSS());
                figureCollection.addFigure(featureFigure);
            }
            featureFigure.forceRegeneration();
        }

        Collection<SimpleFeatureFigure> remainingFigures = figureMap.values();
        figureCollection.removeFigures(remainingFigures.toArray(new Figure[remainingFigures.size()]));

    }

    public SimpleFeatureFigureFactory getFigureFactory() {
        return figureFactory;
    }

    public FigureCollection getFigureCollection() {
        return figureCollection;
    }

    @Override
    protected Rectangle2D getLayerModelBounds() {
        if (figureCollection.getFigureCount() == 0) {
            return null;
        } else {
            return figureCollection.getBounds();
        }
    }

    @Override
    protected void renderLayer(Rendering rendering) {
        final Graphics2D g2d = rendering.getGraphics();
        final Object antiAliasing = g2d.getRenderingHint(RenderingHints.KEY_ANTIALIASING);
        final Object textAntiAliasing = g2d.getRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING);

        g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        g2d.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);

        try {
            figureCollection.draw(rendering);
        } finally {
            g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, antiAliasing);
            g2d.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, textAntiAliasing);
        }
    }

    private class VectorDataChangeHandler extends ProductNodeListenerAdapter {

        @Override
        public void nodeChanged(ProductNodeEvent event) {
            if (event.getSourceNode() == VectorDataLayer.this.vectorDataNode) {
                if (ProductNode.PROPERTY_NAME_NAME.equals(event.getPropertyName())) {
                    setName(VectorDataLayer.this.vectorDataNode.getName());
                }
                if (VectorDataNode.PROPERTY_NAME_FEATURE_COLLECTION.equals(event.getPropertyName())) {
                    if (!adjustingFigures) {
                        updateFigureCollection();
                        fireLayerDataChanged(null); // todo - compute changed modelRegion instead of passing null (nf)
                    }
                }
            }
        }
    }

    private class FigureCollectionChangeListener extends AbstractFigureChangeListener {
        @Override
        public void figureChanged(FigureChangeEvent event) {
            final Figure sourceFigure = event.getSourceFigure();
            if (sourceFigure instanceof SimpleFeatureFigure) {
                SimpleFeatureFigure featureFigure = (SimpleFeatureFigure) sourceFigure;
                final SimpleFeature[] features = {featureFigure.getSimpleFeature()};
                try {
                    adjustingFigures = true;
                    VectorDataLayer.this.vectorDataNode.fireFeatureCollectionChanged(features, features);
                    fireLayerDataChanged(null); // todo - compute changed modelRegion instead of passing null (nf)
                } finally {
                    adjustingFigures = false;
                }
            }
        }

        @Override
        public void figuresAdded(FigureChangeEvent event) {

        }

        @Override
        public void figuresRemoved(FigureChangeEvent event) {

        }
    }
}