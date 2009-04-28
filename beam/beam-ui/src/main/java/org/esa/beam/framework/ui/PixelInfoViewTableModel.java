/*
 * $Id: PixelInfoViewTableModel.java,v 1.1 2009-04-28 14:17:18 lveci Exp $
 *
 * Copyright (C) 2009 by Brockmann Consult (info@brockmann-consult.de)
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
package org.esa.beam.framework.ui;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import javax.swing.table.AbstractTableModel;

/**
 *
 * @author Marco Zuehlke
 * @version $Revision: 1.1 $ $Date: 2009-04-28 14:17:18 $
 * @since BEAM 4.5.2
 */
class PixelInfoViewTableModel extends AbstractTableModel {

    private final String[] columnNames;
    private final List<String> names;
    private final List<String> values;
    private final List<String> units;

    public PixelInfoViewTableModel(String[] columnNames) {
        this.columnNames = columnNames;
        names = Collections.synchronizedList(new ArrayList<String>(32));
        values = Collections.synchronizedList(new ArrayList<String>(32));
        units = Collections.synchronizedList(new ArrayList<String>(32));
    }
    
    @Override
    public String getColumnName(int column) {
        return columnNames[column];
    }
    
    @Override
    public int getColumnCount() {
        return columnNames.length;
    }

    @Override
    public int getRowCount() {
        return names.size();
    }

    @Override
    public Object getValueAt(int rowIndex, int columnIndex) {
        if (columnIndex == 0) {
            return names.get(rowIndex);
        } else if (columnIndex == 1) {
            return values.get(rowIndex);
        } else if (columnIndex == 2) {
            return units.get(rowIndex);
        }
        return "";
    }
    
    public void addRow(String name, String value, String unit) {
        synchronized (this) {
            names.add(name);
            values.add(value);
            units.add(unit);
        }
    }
    
    public void updateValue(String aValue, int row) {
        values.set(row, aValue);
    }
    
    public void clear() {
        synchronized (this) {
            names.clear();
            values.clear();
            units.clear();
        }        
    }
}
