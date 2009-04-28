/*
 * $Id: FileArrayEditor.java,v 1.1 2009-04-28 14:17:18 lveci Exp $
 *
 * Copyright (C) 2002 by Brockmann Consult (info@brockmann-consult.de)
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
package org.esa.beam.framework.ui.io;

import org.esa.beam.framework.ui.GridBagUtils;
import org.esa.beam.framework.ui.UIUtils;
import org.esa.beam.framework.ui.tool.ToolButtonFactory;
import org.esa.beam.util.Guardian;
import org.esa.beam.util.SystemUtils;
import org.esa.beam.util.io.BeamFileChooser;

import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.util.ArrayList;
import java.util.List;

/**
 * An UI-Component which represents a file list with the ability to add and remove files.
 */
public final class FileArrayEditor {

    private static final Dimension _listPreferredSize = new Dimension(500, 200);
    private JPanel _basePanel;
    private JList _listComponent;
    private JFileChooser _fileDialog;
    private final List _fileList;
    private FileArrayEditorListener _listener;
    private final EditorParent _parent;
    private String _label;

    /**
     * Constructs the object with default values
     *
     * @param parent the parent editor
     * @param label  the label for this editor
     */
    public FileArrayEditor(final EditorParent parent, String label) {
        Guardian.assertNotNullOrEmpty("label", label);
        _parent = parent;
        _label = label;
        _basePanel = null;
        _listener = null;
        _fileList = new ArrayList();
    }

    /**
     * Retrieves the UI of the editor. Returns the base component
     */
    public JComponent getUI() {
        if (_basePanel == null) {
            createUi();
        }

        return _basePanel;
    }

    /**
     * Sets the list of files to be edited. The list currently held is overwritten.
     *
     * @param files <code>List</code> of <code>File</code>s to be set
     */
    public void setFiles(final List files) {
        Guardian.assertNotNull("files", files);
        _fileList.clear();
        _fileList.addAll(files);
        _listComponent.setListData(_fileList.toArray());
        notifyListener();
    }

    /**
     * Retrieves the list of files currently edited
     *
     * @return a <code>List</code> of currently edited <code>File</code>s
     */
    public List getFiles() {
        return _fileList;
    }

    /**
     * Sets the listener for this class
     *
     * @param listener the listener to associate with this editor
     */
    public void setListener(final FileArrayEditorListener listener) {
        _listener = listener;
    }

    ///////////////////////////////////////////////////////////////////////////
    ////// END OF PUBLIC
    ///////////////////////////////////////////////////////////////////////////


    /**
     * Creates the user interface
     */
    private void createUi() {
        // the label
        final JLabel label = new JLabel(_label + ":");
        setName(label, _label);

        // the list
        _listComponent = new JList();
        setName(_listComponent, _label);
        JScrollPane scrollPane = new JScrollPane(_listComponent);
        setName(scrollPane, _label);
        scrollPane.setPreferredSize(_listPreferredSize);

        // the add button
        final JButton addButton = (JButton) ToolButtonFactory.createButton(
                UIUtils.loadImageIcon("icons/Plus16.gif"),
                false);
        setName(addButton, "addButton");
        addButton.addActionListener(
                new ActionListener() {

                    public void actionPerformed(final ActionEvent e) {
                        onAddButton();
                    }
                });
        // the remove button
        final JButton removeButton = (JButton) ToolButtonFactory.createButton(
                UIUtils.loadImageIcon("icons/Minus16.gif"), false);
        setName(removeButton, "removeButton");
        removeButton.addActionListener(
                new ActionListener() {

                    public void actionPerformed(final ActionEvent e) {
                        onRemoveButton();
                    }
                });

        // the button panel
        final JPanel buttonPanel = new JPanel();
        setName(buttonPanel, _label);
        buttonPanel.add(addButton);
        buttonPanel.add(removeButton);

        // the base panel
        _basePanel = GridBagUtils.createDefaultEmptyBorderPanel();
        setName(_basePanel, _label);
        final GridBagConstraints gbc = GridBagUtils.createConstraints(null);
        gbc.anchor = GridBagConstraints.WEST;
        gbc.weightx = 1;

        gbc.gridy++;
        _basePanel.add(label, gbc);
        gbc.anchor = GridBagConstraints.EAST;
        _basePanel.add(buttonPanel, gbc);

        gbc.gridy++;
        gbc.anchor = GridBagConstraints.WEST;
        gbc.gridwidth = 2;
        gbc.fill = GridBagConstraints.BOTH;
        gbc.weightx = 1;
        gbc.weighty = 1;
        _basePanel.add(scrollPane, gbc);

    }

    private void setName(final Component comp, String name) {
        comp.setName(name);
    }

    /**
     * Callback invoked by the add button
     */
    private void onAddButton() {
        _fileDialog = getFileDialogSafe();
        final File userInputDir = _parent.getUserInputDir();
        final int retVal;

        _fileDialog.setCurrentDirectory(userInputDir);
        retVal = _fileDialog.showOpenDialog(_basePanel);

        if (retVal == JFileChooser.APPROVE_OPTION) {
            File[] selected = _fileDialog.getSelectedFiles();

            for (int n = 0; n < selected.length; n++) {
                _fileList.add(selected[n]);
            }

            _listComponent.setListData(_fileList.toArray());
            notifyListener();
            _parent.setUserInputDir(_fileDialog.getCurrentDirectory());
        }
    }

    /**
     * Callback invoked by the remove button
     */
    private void onRemoveButton() {
        final Object[] toRemove = _listComponent.getSelectedValues();
        for (int n = 0; n < toRemove.length; n++) {
            _fileList.remove(toRemove[n]);
        }
        _listComponent.setListData(_fileList.toArray());
        notifyListener();
    }

    /**
     * Retrieves the file chooser object. If none is present, an object is constructed
     */
    private JFileChooser getFileDialogSafe() {
        if (_fileDialog == null) {
            _fileDialog = new BeamFileChooser();
            _fileDialog.setFileSelectionMode(JFileChooser.FILES_ONLY);
            _fileDialog.setCurrentDirectory(SystemUtils.getUserHomeDir());
            _fileDialog.setMultiSelectionEnabled(true);
        }

        return _fileDialog;
    }

    /**
     * Calls the listener about changes - if necessary
     */
    private void notifyListener() {
        if ((_listener != null)) {
            _listener.updatedList((File[]) _fileList.toArray(new File[0]));
        }
    }

    public static interface EditorParent {

        File getUserInputDir();

        void setUserInputDir(File newDir);
    }

    public interface FileArrayEditorListener {

        public void updatedList(File[] files);
    }
}
