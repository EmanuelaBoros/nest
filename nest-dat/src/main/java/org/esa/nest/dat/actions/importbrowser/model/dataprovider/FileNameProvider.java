/*
 * Created on: 	19.09.2005
 * Created by:	Marco Peters
 * File: 		FileNameProvider.java 
 */
package org.esa.nest.dat.actions.importbrowser.model.dataprovider;

import java.awt.*;
import java.io.IOException;
import java.util.Comparator;

import javax.swing.JLabel;
import javax.swing.JTable;
import javax.swing.table.DefaultTableCellRenderer;
import javax.swing.table.TableColumn;

import org.esa.nest.dat.actions.importbrowser.model.Repository;
import org.esa.nest.dat.actions.importbrowser.model.RepositoryEntry;

public class FileNameProvider implements DataProvider {

    private final Comparator _fileNameComparator = new FileNameComparator();
    private TableColumn _fileNameColumn;

    public boolean mustCreateData(final RepositoryEntry entry, final Repository repository) {
        return false;
    }

    public void createData(final RepositoryEntry entry, final Repository repository)
            throws IOException {

    }

    public Object getData(final RepositoryEntry entry, final Repository repository)
            throws IOException {
        return entry.getProductFile().getName();
    }

    public Comparator getComparator() {
        return _fileNameComparator;
    }

    public void cleanUp(final RepositoryEntry entry, final Repository repository) {
    }

    public TableColumn getTableColumn() {
        if(_fileNameColumn == null) {
            _fileNameColumn = new TableColumn();
            _fileNameColumn.setHeaderValue("File Name");
            _fileNameColumn.setPreferredWidth(150);
            _fileNameColumn.setResizable(true);
            _fileNameColumn.setCellRenderer(new FileNameCellRenderer());
        }
        return _fileNameColumn;
    }

    private static class FileNameCellRenderer extends DefaultTableCellRenderer {

        @Override
        public Component getTableCellRendererComponent(final JTable table,
                                                       final Object value,
                                                       final boolean isSelected,
                                                       final boolean hasFocus,
                                                       final int row, final int column) {
            final JLabel jlabel = (JLabel) super
                    .getTableCellRendererComponent(table, value, isSelected, hasFocus,
                                                   row, column);

            jlabel.setFont(jlabel.getFont().deriveFont(Font.BOLD));
            String valStr = (String) value;
            jlabel.setToolTipText(valStr);
            return jlabel;
        }
    }

    private static class FileNameComparator implements Comparator {

        public int compare(final Object o1, final Object o2) {
            if(o1 == o2) {
                return 0;
            }
            if (o1 == null) {
                return -1;
            } else if(o2 == null) {
                return 1;
            }
                
            final String s1 = (String) o1;
            final String s2 = (String) o2;

            return s1.compareTo(s2);
        }
    }
}
