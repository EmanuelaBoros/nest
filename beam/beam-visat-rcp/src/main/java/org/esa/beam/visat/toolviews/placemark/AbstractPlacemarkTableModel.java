package org.esa.beam.visat.toolviews.placemark;

import com.bc.ceres.core.ProgressMonitor;
import org.esa.beam.framework.datamodel.Band;
import org.esa.beam.framework.datamodel.GeoPos;
import org.esa.beam.framework.datamodel.Pin;
import org.esa.beam.framework.datamodel.PixelPos;
import org.esa.beam.framework.datamodel.PlacemarkDescriptor;
import org.esa.beam.framework.datamodel.Product;
import org.esa.beam.framework.datamodel.ProductNodeEvent;
import org.esa.beam.framework.datamodel.ProductNodeListenerAdapter;
import org.esa.beam.framework.datamodel.TiePointGrid;
import org.esa.beam.util.math.MathUtils;

import javax.swing.table.DefaultTableModel;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;


public abstract class AbstractPlacemarkTableModel extends DefaultTableModel {

    private final PlacemarkDescriptor placemarkDescriptor;

    private Product product;
    private Band[] selectedBands;
    private TiePointGrid[] selectedGrids;

    private final PlacemarkListener placemarkListener;
    private final ArrayList<Pin> placemarkList;

    public AbstractPlacemarkTableModel(PlacemarkDescriptor placemarkDescriptor, Product product, Band[] selectedBands,
                                       TiePointGrid[] selectedGrids) {
        this.placemarkDescriptor = placemarkDescriptor;
        this.product = product;
        initSelectedBands(selectedBands);
        initSelectedGrids(selectedGrids);
        placemarkList = new ArrayList<Pin>(10);
        placemarkListener = new PlacemarkListener();
        if (product != null) {
            product.addProductNodeListener(placemarkListener);
        }
        initPlacemarkList(product);
    }

    public Pin[] getPlacemarks() {
        return placemarkList.toArray(new Pin[placemarkList.size()]);
    }

    public Pin getPlacemarkAt(int modelRow) {
        return placemarkList.get(modelRow);
    }

    public PlacemarkDescriptor getPlacemarkDescriptor() {
        return placemarkDescriptor;
    }

    public Product getProduct() {
        return product;
    }

    public void setProduct(Product product) {
        if (this.product == product) {
            return;
        }
        if (this.product != null) {
            this.product.removeProductNodeListener(placemarkListener);
        }
        this.product = product;
        if (this.product != null) {
            this.product.addProductNodeListener(placemarkListener);
        }

        placemarkList.clear();
        initPlacemarkList(this.product);
        selectedBands = null;
        selectedGrids = null;
        fireTableStructureChanged();
    }

    public Band[] getSelectedBands() {
        return selectedBands;
    }

    public void setSelectedBands(Band[] selectedBands) {
        this.selectedBands = selectedBands != null ? selectedBands : new Band[0];
        fireTableStructureChanged();
    }

    public TiePointGrid[] getSelectedGrids() {
        return selectedGrids;
    }

    public void setSelectedGrids(TiePointGrid[] selectedGrids) {
        this.selectedGrids = selectedGrids != null ? selectedGrids : new TiePointGrid[0];
        fireTableStructureChanged();
    }

    public boolean addPlacemark(Pin placemark) {
        if (getProduct() != null && placemarkList.add(placemark)) {
            fireTableDataChanged();
            return true;
        }
        return false;
    }

    public boolean removePlacemark(Pin placemark) {
        if (getProduct() != null && placemarkList.remove(placemark)) {
            fireTableDataChanged();
            return true;
        }
        return false;
    }

    public void removePlacemarkAt(int index) {
        placemarkList.remove(index);
    }

    public abstract String[] getStandardColumnNames();

    @Override
    public abstract boolean isCellEditable(int rowIndex, int columnIndex);

    protected abstract Object getStandardColumnValueAt(int rowIndex, int columnIndex);

    @Override
    public int getRowCount() {
        if (placemarkList == null) {
            return 0;
        }
        return placemarkList.size();
    }

    @Override
    public int getColumnCount() {
        int count = getStandardColumnNames().length;
        if (selectedBands != null) {
            count += selectedBands.length;
        }
        if (selectedGrids != null) {
            count += selectedGrids.length;
        }
        return count;
    }

    @Override
    public String getColumnName(int columnIndex) {
        if (columnIndex < getStandardColumnNames().length) {
            return getStandardColumnNames()[columnIndex];
        }
        int newIndex = columnIndex - getStandardColumnNames().length;
        if (newIndex < getNumSelectedBands()) {
            return selectedBands[newIndex].getName();
        }
        newIndex -= getNumSelectedBands();
        if (selectedGrids != null && newIndex < selectedGrids.length) {
            return selectedGrids[newIndex].getName();
        }
        return "?";
    }

    @Override
    public Class getColumnClass(int columnIndex) {
        if (columnIndex >= 0 && columnIndex < getStandardColumnNames().length - 1) {
            return Float.class;
        }
        return Object.class;
    }

    @Override
    public Object getValueAt(int rowIndex, int columnIndex) {
        if (product != null) {
            if (columnIndex < getStandardColumnNames().length) {
                return getStandardColumnValueAt(rowIndex, columnIndex);
            } else {
                final Pin placemark = placemarkList.get(rowIndex);
                int index = columnIndex - getStandardColumnNames().length;
                PixelPos pixelPos = placemark.getPixelPos();
                if (pixelPos == null) {
                    return "No-data";
                }

                final int x = MathUtils.floorInt(pixelPos.getX());
                final int y = MathUtils.floorInt(pixelPos.getY());
                final int width = product.getSceneRasterWidth();
                final int height = product.getSceneRasterHeight();

                if (x < 0 || x >= width || y < 0 || y >= height) {
                    return "No-data";
                }

                if (index < getNumSelectedBands()) {
                    final Band band = selectedBands[index];
                    if (band.isPixelValid(x, y)) {
                        float[] value = null;
                        try {
                            value = band.readPixels(x, y, 1, 1, value, ProgressMonitor.NULL);
                            return value[0];
                        } catch (IOException e) {
                            return "I/O-error";
                        }
                    } else {
                        return "No-data";
                    }
                }
                index -= getNumSelectedBands();
                if (index < selectedGrids.length) {
                    final TiePointGrid grid = selectedGrids[index];
                    float[] value = null;
                    try {
                        value = grid.readPixels(x, y, 1, 1, value, ProgressMonitor.NULL);
                        return value[0];
                    } catch (IOException e) {
                        return "I/O-error";
                    }
                }
            }
        }

        return "";
    }

    @Override
    public void setValueAt(Object value, int rowIndex, int columnIndex) {
        if (value == null) {
            return;
        }
        if (columnIndex < getStandardColumnNames().length) {
            Pin placemark = placemarkList.get(rowIndex);
            if (columnIndex == 0) {
                if (value instanceof Float) {
                    float pixelY;
                    if (placemark.getPixelPos() == null) {
                        pixelY = -1;
                    } else {
                        pixelY = placemark.getPixelPos().y;
                    }
                    placemark.setPixelPos(new PixelPos((Float) value, pixelY));
                    GeoPos geoPos = placemarkDescriptor.updateGeoPos(product.getGeoCoding(),
                                                                     placemark.getPixelPos(),
                                                                     placemark.getGeoPos());
                    placemark.setGeoPos(geoPos);
                }
            } else if (columnIndex == 1) {
                if (value instanceof Float) {
                    float pixelX;
                    if (placemark.getPixelPos() == null) {
                        pixelX = -1;
                    } else {
                        pixelX = placemark.getPixelPos().x;
                    }
                    placemark.setPixelPos(new PixelPos(pixelX, (Float) value));
                    GeoPos geoPos = placemarkDescriptor.updateGeoPos(product.getGeoCoding(),
                                                                     placemark.getPixelPos(),
                                                                     placemark.getGeoPos());
                    placemark.setGeoPos(geoPos);
                }
            } else if (columnIndex == 2) {
                if (value instanceof Float) {
                    float lat;
                    if (placemark.getGeoPos() == null) {
                        lat = Float.NaN;
                    } else {
                        lat = placemark.getGeoPos().lat;
                    }
                    placemark.setGeoPos(new GeoPos(lat, (Float) value));
                    PixelPos pixelPos = placemarkDescriptor.updatePixelPos(product.getGeoCoding(),
                                                                           placemark.getGeoPos(),
                                                                           placemark.getPixelPos());
                    placemark.setPixelPos(pixelPos);
                }
            } else if (columnIndex == 3) {
                if (value instanceof Float) {
                    float lon;
                    if (placemark.getGeoPos() == null) {
                        lon = Float.NaN;
                    } else {
                        lon = placemark.getGeoPos().lon;
                    }
                    placemark.setGeoPos(new GeoPos((Float) value, lon));
                    PixelPos pixelPos = placemarkDescriptor.updatePixelPos(product.getGeoCoding(),
                                                                           placemark.getGeoPos(),
                                                                           placemark.getPixelPos());
                    placemark.setPixelPos(pixelPos);
                }
            } else if (columnIndex == getStandardColumnNames().length - 1) {
                String strValue = value.toString();
                placemark.setLabel(strValue);
            } else {
                throw new IllegalStateException(
                        "Column[" + columnIndex + "] '" + getColumnName(columnIndex) + "' is not editable");
            }
        }
    }

    public void dispose() {
        if (product != null) {
            product.removeProductNodeListener(placemarkListener);
        }
        selectedBands = null;
        selectedGrids = null;
        placemarkList.clear();
    }

    private void initSelectedBands(Band[] selectedBands) {
        this.selectedBands = selectedBands != null ? selectedBands : new Band[0];
    }

    private void initSelectedGrids(TiePointGrid[] selectedGrids) {
        this.selectedGrids = selectedGrids != null ? selectedGrids : new TiePointGrid[0];
    }

    private void initPlacemarkList(Product product) {
        if (product != null) {
            Pin[] placemarks = placemarkDescriptor.getPlacemarkGroup(product).toArray(new Pin[0]);
            placemarkList.addAll(Arrays.asList(placemarks));
        }
    }

    private int getNumSelectedBands() {
        return selectedBands != null ? selectedBands.length : 0;
    }

    private class PlacemarkListener extends ProductNodeListenerAdapter {

        @Override
        public void nodeChanged(ProductNodeEvent event) {
            fireTableDataChanged(event);
        }

        private void fireTableDataChanged(ProductNodeEvent event) {
            if (event.getSourceNode() instanceof Pin) {
                Pin placemark = (Pin) event.getSourceNode();
                if (placemarkList.contains(placemark)) {
                    AbstractPlacemarkTableModel.this.fireTableDataChanged();
                }
            }
        }
    }
}
