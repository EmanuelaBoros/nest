/*
 * $Id: BeamFileFilter.java,v 1.2 2009-07-07 00:27:41 lveci Exp $
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
package org.esa.beam.util.io;

import org.esa.beam.util.StringUtils;

import javax.swing.filechooser.FileFilter;
import javax.swing.JFileChooser;
import java.io.File;
import java.util.ArrayList;

/**
 * A <code>FileFilter</code> with file extensions support.
 *
 * @author Norman Fomferra
 * @version $Revision: 1.2 $  $Date: 2009-07-07 00:27:41 $
 */
public class BeamFileFilter extends FileFilter {

    private String formatName;
    private String[] extensions;
    private String description;

    public BeamFileFilter() {
    }

    public BeamFileFilter(String formatName, String extension, String description) {
        this(formatName, StringUtils.toStringArray(extension, ","), description);
    }

    public BeamFileFilter(String formatName, String[] extensions, String description) {
        setFormatName(formatName);
        setExtensions(extensions);
        setDescription(description);
    }

    public String getFormatName() {
        return formatName;
    }

    public void setFormatName(String formatName) {
        this.formatName = formatName;
    }

    /**
     * Returns whether or not this file filter has extensions.
     *
     * @return <code>true</code> if so
     */
    public boolean hasExtensions() {
        return extensions != null && extensions.length > 0;
    }

    /**
     * Returns the default extension. The default extension is the first entry in the array returned by the
     * <code>getExtensions</code> method.
     *
     * @return the default extension or <code>null</code> if no extensions have bees specified.
     *
     * @see #getExtensions
     */
    public String getDefaultExtension() {
        return hasExtensions() ? getExtensions()[0] : null;
    }

    /**
     * Returns the accepted extensions of this filter. For example: <code>{".jpg", ".gif", ".png"}</code>.
     *
     * @return The array of extensions.
     *
     * @see #setExtensions
     */
    public String[] getExtensions() {
        return extensions;
    }

    /**
     * Sets the accepted extensions of this filter. For example: <code>{".jpg", ".gif", ".png"}</code>.
     *
     * @param extensions The array of extensions.
     *
     * @see #getExtensions
     */
    public void setExtensions(String[] extensions) {
        if (extensions != null) {
            ArrayList<String> extensionList = new ArrayList<String>();
            for (final String extension : extensions) {
                if (extension.startsWith(".")) {
                    extensionList.add(extension);
                } else if (extension.trim().length() > 0) {
                    extensionList.add("." + extension);
                }
            }
            this.extensions = extensionList.toArray(new String[extensionList.size()]);
        } else {
            this.extensions = null;
        }
    }

    /**
     * Returns the description of this filter. For example: <code>"JPEG Images (*.jpg,*.jpeg)"</code>.
     *
     * @see javax.swing.filechooser.FileView#getTypeDescription(java.io.File)
     */
    @Override
    public String getDescription() {
        return description;
    }

    /**
     * Returns the description of this filter. For example: <code>"JPEG Images (*.jpg,*.jpeg)"</code>. If the extension
     * list is missing in the description text, it is automatically appended.
     *
     * @param description The description, must not be null.
     *
     * @see #getDescription
     */
    public void setDescription(String description) {
        if (hasExtensions() && !description.endsWith(")")) {
            StringBuffer sb = new StringBuffer(description);
            sb.append(" (");
            for (int i = 0; i < extensions.length; i++) {
                if (i > 0) {
                    sb.append(",");
                }
                sb.append("*");
                if (extensions[i] != null) {
                    sb.append(extensions[i]);
                }
            }
            sb.append(")");
            this.description = sb.toString();
        } else {
            this.description = description;
        }
    }

    /**
     * Utility method which checks the extension the given file.
     *
     * @param file the file
     *
     * @return <code>true</code> if the given file path ends with one of the registered extensions, <code>false</code>
     *         otherwise.
     */
    public boolean checkExtension(File file) {
        return file != null && checkExtension(file.getName());
    }

    /**
     * Utility method which checks the extension the given filename.
     *
     * @param filename the file name
     *
     * @return <code>true</code> if the given file name ends with one of the registered extensions, <code>false</code>
     *         otherwise.
     */
    public boolean checkExtension(String filename) {
        if (filename != null) {
            filename = filename.toLowerCase();
            for (String extension : extensions) {
                extension = extension.toLowerCase();
                if (filename.endsWith(extension)) {
                    return true;
                }
            }
        }
        return false;
    }

    /**
     * Tests whether or not the given file is accepted by this filter. The default implementation returns
     * <code>true</code> if the given file is a directory or the path string ends with one of the registered extensions.
     * if no extension are defined, the method always returns <code>true</code>
     *
     * @param file the file to be or not be accepted.
     *
     * @return <code>true</code> if given file is accepted by this filter
     */
    @Override
    public boolean accept(File file) {
        if (!hasExtensions()) {
            return true;
        }

        // directories are accepted right away
        if (file.isDirectory()) {
            return true;
        }

        // otherwise name must end with one of the extensions
        return checkExtension(file);
    }

    /**
     * Checks if the given directory represents a compound document.
     * If so, we don't want the user to descend into it when using the
     * {@link org.esa.beam.util.io.BeamFileChooser}.
     * The default implementation returns {@code false}.
     * Clients may override.
     *
     * @param dir The directory to check.
     *
     * @return {@code true} If the given directory represents a compound document.
     * @since BEAM 4.6.1
     */
    public boolean isCompoundDocument(File dir) {
        return false;
    }

    /**
     * Gets the file selection mode for the {@link org.esa.beam.util.io.BeamFileChooser} if this filter is used.
     * The default implementation returns {@link FileSelectionMode#FILES_ONLY}.
     * Clients may override.
     *
     * @return {@code true} if the user can also select directories using this filter.
     * @since BEAM 4.6.1
     */
    public FileSelectionMode getFileSelectionMode() {
        return FileSelectionMode.FILES_ONLY;
    }

    /**
     * File selection modes.
     */
    public enum FileSelectionMode {
        /** Instruction to display only files. */
        FILES_ONLY(JFileChooser.FILES_ONLY),

        /** Instruction to display only directories. */
        DIRECTORIES_ONLY(JFileChooser.DIRECTORIES_ONLY),

        /** Instruction to display both files and directories. */
        FILES_AND_DIRECTORIES (JFileChooser.FILES_AND_DIRECTORIES);

        private final int value;

        FileSelectionMode(int value) {
            this.value = value;
        }

        public int getValue() {
            return value;
        }
    }
}
