/*
 * Created on: 	19.09.2005
 * Created by:	Marco Peters
 * File: 		ProductSizeProvider.java 
 */
package org.esa.nest.dat.actions.importbrowser.model.dataprovider;

import org.esa.beam.framework.datamodel.Product;
import org.esa.nest.dat.actions.importbrowser.model.Repository;
import org.esa.nest.dat.actions.importbrowser.model.RepositoryEntry;
import org.esa.nest.dataio.ReaderUtils;

import javax.swing.*;
import javax.swing.table.DefaultTableCellRenderer;
import javax.swing.table.TableColumn;
import java.awt.*;
import java.io.IOException;
import java.util.Comparator;

public class ProductSizeProvider implements DataProvider {

    private TableColumn _fileSizeColumn;
    private final Comparator _productSizeComparator = new ProductSizeComparator();

    public boolean mustCreateData(final RepositoryEntry entry,final  Repository repository) {
        return false;
    }

    public void createData(final RepositoryEntry entry,final  Repository repository)
            throws IOException {
    }

    public Object getData(final RepositoryEntry entry,final  Repository repository)
            throws IOException {
        float prodSize = entry.getProductSize();
        Product product = entry.getProduct();
        if(prodSize < 1 && product != null)
            prodSize = ReaderUtils.getTotalSize(product);
        return prodSize;
    }

    public Comparator getComparator() {
        return _productSizeComparator;
    }

    public void cleanUp(final RepositoryEntry entry,final  Repository repository) {
    }

    public TableColumn getTableColumn() {
        if (_fileSizeColumn == null) {
            _fileSizeColumn = new TableColumn();
            _fileSizeColumn.setHeaderValue("Size");
            _fileSizeColumn.setPreferredWidth(70);
            _fileSizeColumn.setResizable(true);
            _fileSizeColumn.setCellRenderer(new FileSizeCellRenderer());
        }
        return _fileSizeColumn;
    }

    private static class FileSizeCellRenderer extends DefaultTableCellRenderer {

        public Component getTableCellRendererComponent(final JTable table,
                                                       final Object value,
                                                       final boolean isSelected,
                                                       final boolean hasFocus,
                                                       final int row, final int column) {
            final JLabel jlabel = (JLabel) super
                    .getTableCellRendererComponent(table, value, isSelected, hasFocus,
                                                   row, column);

            jlabel.setFont(jlabel.getFont().deriveFont(Font.BOLD));
            jlabel.setHorizontalAlignment(JLabel.CENTER);
            setText(String.format("%1$.2f MB", new Object[] {value }));
            return jlabel;
        }
    }

    private static class ProductSizeComparator implements Comparator {

        public int compare(final Object o1, final Object o2) {
            if(o1 == o2) {
                return 0;
            }
            if (o1 == null) {
                return -1;
            } else if(o2 == null) {
                return 1;
            }

            final Float f1 = (Float) o1;
            final Float f2 = (Float) o2;

            return f1.compareTo(f2);
        }
    }
}
