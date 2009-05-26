/*
 * $Id: SystemUtils.java,v 1.2 2009-05-26 19:44:50 lveci Exp $
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
package org.esa.beam.util;

import com.bc.ceres.core.runtime.internal.RuntimeActivator;
import org.esa.beam.util.logging.BeamLogManager;

import javax.swing.UIManager;
import java.awt.Image;
import java.awt.Toolkit;
import java.awt.datatransfer.Clipboard;
import java.awt.datatransfer.DataFlavor;
import java.awt.datatransfer.StringSelection;
import java.awt.datatransfer.Transferable;
import java.awt.datatransfer.UnsupportedFlavorException;
import java.io.File;
import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.net.URI;
import java.net.URISyntaxException;
import java.net.URL;
import java.net.URLDecoder;
import java.security.CodeSource;
import java.text.MessageFormat;
import java.util.NoSuchElementException;
import java.util.ServiceLoader;
import java.util.StringTokenizer;

/**
 * A collection of (BEAM-) system level functions.
 * <p/>
 * <p> All functions have been implemented with extreme caution in order to provide a maximum performance.
 *
 * @author Norman Fomferra
 * @author Sabine Embacher
 * @version $Revision: 1.2 $ $Date: 2009-05-26 19:44:50 $
 */
public class SystemUtils {

    /**
     * The URL string to the BEAM home page.
     */
    public static final String BEAM_HOME_PAGE = "http://envisat.esa.int/beam/";

    public static final String BEAM_HOME_PROPERTY_NAME = "beam.home";
    public static final String LAX_INSTALL_DIR_PROPERTY_NAME = "lax.root.install.dir";
    /**
     * SYSTEM_DEPENDENT_LINE_SEPARATOR
     */
    public static final String LS = System.getProperty("line.separator");

    private static final char _URL_DIR_SEPARATOR_CHAR = '/';
    public static final int LL_DEBUG = 10;
    public static final int LL_INFO = 20;
    public static final int LL_WARNING = 30;
    public static final int LL_ERROR = 40;
    public static final String LLS_DEBUG = "DEBUG";
    public static final String LLS_INFO = "INFO";
    public static final String LLS_WARNING = "WARNING";
    public static final String LLS_ERROR = "ERROR";

    /**
     * Name of BEAM's extensions directory.
     */
    public static final String EXTENSION_DIR_NAME = "extensions";

    /**
     * Name of BEAM's auxdata directory.
     */
    public static final String AUXDATA_DIR_NAME = "auxdata";
    public static final String CACHE_DIR_NAME = "cache";
    public static final String BEAM_PLUGIN_PATH_PROPERTY_NAME = "beam.plugin.path";

    private static final String _H5_CLASS_NAME = "ncsa.hdf.hdf5lib.H5";
    private static final String _H4_CLASS_NAME = "ncsa.hdf.hdflib.HDFLibrary";
    private static final String FILE_PROTOCOLL_PREFIX = "file:";
    private static final String JAR_PROTOCOL_PREFIX = "jar:";

    /**
     * Gets the current user's name, or the string <code>"unknown"</code> if the the user's name cannot be determined.
     *
     * @return the current user's name, never <code>null</code>
     */
    public static String getUserName() {
        return System.getProperty("user.name", "unknown");
    }

    /**
     * Gets the current user's home directory, or the directory pointed to by '.' if the user's actual home directory
     * cannot be determined.
     *
     * @return the current working directory, never <code>null</code>
     */
    public static File getUserHomeDir() {
        return new File(System.getProperty("user.home", "."));
    }

    /**
     * Gets the current user's application data directory.
     *
     * @return the current user's application data directory
     *
     * @since BEAM 4.2
     */
    public static File getApplicationDataDir() {
        return getApplicationDataDir(false);
    }

    /**
     * Optionally creates and returns the current user's application data directory.
     *
     * @param force if true, the directory will be created if it didn't exist before
     *
     * @return the current user's application data directory
     *
     * @since BEAM 4.2
     */
    public static File getApplicationDataDir(boolean force) {
        String contextId = getApplicationContextId();
        final File dir = new File(getUserHomeDir(), "." + contextId);
        if (force && !dir.exists()) {
            dir.mkdirs();
        }
        return dir;
    }

    private static String getApplicationContextId() {
        String contextId = null;
        if (RuntimeActivator.getInstance() != null
                && RuntimeActivator.getInstance().getModuleContext() != null) {
            contextId = RuntimeActivator.getInstance().getModuleContext().getRuntimeConfig().getContextId();
        }
        if (contextId == null) {
            contextId = System.getProperty("ceres.context", "beam");
        }
        return contextId;
    }

    /**
     * Gets the current working directory, or the directory pointed to by '.' if the actual working directory cannot be
     * determined.
     *
     * @return the current working directory, never <code>null</code>
     */
    public static File getCurrentWorkingDir() {
        return new File(System.getProperty("user.dir", "."));
    }

    /**
     * Gets all files (class directory & JAR file pathes) given in the current class path of the Java runtime which
     * loaded this class.
     * <p/>
     * <p> The files pathes returned are either relative or absolute, just as they where defined for the runtime's class
     * path.
     *
     * @return all files in the current class path, never <code>null</code>
     */
    public static File[] getClassPathFiles() {
        String classPath = System.getProperty("java.class.path");
        if (classPath == null) {
            return new File[0];
        }
        StringTokenizer st = new StringTokenizer(classPath, File.pathSeparator);
        File[] files = new File[st.countTokens()];
        try {
            for (int i = 0; i < files.length; i++) {
                files[i] = new File(st.nextToken());
            }
        } catch (NoSuchElementException e) {
            // ignore
        }
        return files;
    }

    /**
     * Gets an application's home directory. The method determines the home directory by retrieving the URL of this
     * class using the method {@link #getApplicationHomeDir(java.net.URL)}.
     *
     * @return an assumption of an application's home directory, never <code>null</code>
     */
    public static File getApplicationHomeDir() {
        final String id = getApplicationContextId();
        final String homeDirPath = System.getProperty(id + ".home");
        if (homeDirPath != null) {
            return new File(homeDirPath);
        }
        // Use fallback
        final URL url = SystemUtils.class.getResource(getClassFileName(SystemUtils.class));
        return getApplicationHomeDir(url);
    }

    /**
     * Extracts an application's home directory from the given URL.
     * <p/>
     * The URL is than scanned for the last occurence of the string <code>&quot;/modules/&quot;</code>.
     * If this succeeds the method returns the absolute
     * (parent) path to the directory which contains <code>modules</code>, which is
     * then assumed to be the requested home directory.
     *
     * @param url the URL
     *
     * @return an assumption of an application's home directory, never <code>null</code>
     *
     * @throws IllegalArgumentException if the given url is <code>null</code>.
     */
    public static File getApplicationHomeDir(final URL url) {
        Guardian.assertNotNull("url", url);
        String path = url.getPath();
        try {
            path = URLDecoder.decode(path, "UTF-8");
        } catch (UnsupportedEncodingException e) {
            // ignored
        }
        path = stripUrlProtocolPrefixes(path);
        path = path.replace(File.separatorChar, '/');
        path = stripClassLibraryPaths(path);
        path = path.replace('/', File.separatorChar);
        return new File(path);
    }

    private static String stripClassLibraryPaths(String path) {
        int pos = path.lastIndexOf("/modules/");
        if (pos >= 0) {
            path = path.substring(0, pos);
        }
        return path;
    }


    /**
     * Retrieves the file name of a class. For example, the string <code>"Date.class"</code> is returned for the
     * class <code>java.util.Date</code>.
     *
     * @param aClass The class.
     *
     * @return the file name of the given class
     *
     * @throws IllegalArgumentException if the given parameter is <code>null</code>.
     */
    public static String getClassFileName(final Class aClass) {
        Guardian.assertNotNull("aClass", aClass);
        final String qualClassName = aClass.getName();
        final int pos = qualClassName.lastIndexOf('.');
        final String className;
        if (pos > 0) {
            className = qualClassName.substring(pos + 1);
        } else {
            className = qualClassName;
        }
        return className + ".class";
    }


    private static String stripUrlProtocolPrefixes(String path) {
        while (true) {
            if (path.startsWith(FILE_PROTOCOLL_PREFIX)) {
                path = path.substring(FILE_PROTOCOLL_PREFIX.length());
            } else if (path.startsWith(JAR_PROTOCOL_PREFIX)) {
                path = path.substring(JAR_PROTOCOL_PREFIX.length());
            } else {
                break;
            }
        }
        return path;
    }

    /**
     * Gets the BEAM Java home directory. The method evaluates the system property <code>org.esa.beam.home</code>. If it
     * is given, it is returned, otherwise <code>getApplicationHomeDir()</code> is returned.
     *
     * @return the BEAM home directory
     */
    public static File getBeamHomeDir() {

        String homeDir = System.getProperty(BEAM_HOME_PROPERTY_NAME);
        if (homeDir != null && homeDir.length() > 0) {
            return new File(homeDir);
        }
        homeDir = System.getProperty(LAX_INSTALL_DIR_PROPERTY_NAME);
        if (homeDir != null && homeDir.length() > 0) {
            return new File(homeDir);
        }

        final URL url = SystemUtils.class.getResource(getClassFileName(SystemUtils.class));
        String path = url.getPath();
        try {
            path = URLDecoder.decode(path, "UTF-8");
        } catch (UnsupportedEncodingException e) {
            // ignored
        }
        path = path.replace(File.separatorChar, '/');
        String beam4Key = "/beam4/";
        int beam4Index = path.indexOf(beam4Key);
        if (beam4Index != -1) {
            path = path.substring(0, beam4Index + beam4Key.length() - 1);
            path = path.replace('/', File.separatorChar);
            return new File(path);
        } else {
            return new File(".").getAbsoluteFile();
        }
    }


    /**
     * Gets the default BEAM cache directory. This is the directory
     * where BEAM stores temporary data.
     *
     * @return the default cache directory
     */
    public static File getDefaultBeamCacheDir() {
        return new File(getApplicationDataDir(), CACHE_DIR_NAME);
    }

    /**
     * Gets the BEAM auxdata directory.
     * This is the directory where auxiliary data is searched.
     * <p>Its value is <code><i>$BEAM_HOME</i>/auxdata</code>.</p>
     *
     * @return the auxdata directory
     *
     * @deprecated in 4.0, use {@link ResourceScanner} instead
     */
    @Deprecated
    public static File getBeamAuxdataDir() {
        CodeSource cs = SystemUtils.class.getProtectionDomain().getCodeSource();
        if (cs != null) {
            URL location = cs.getLocation();
            try {
                URI uri = location.toURI();
                return new File(new File(uri), AUXDATA_DIR_NAME);
            } catch (URISyntaxException e) {
                // ok
            }
        }
        return new File(getBeamHomeDir(), AUXDATA_DIR_NAME);
    }

    /**
     * Replace the separator character '/' with the system-dependent path-separator character.
     *
     * @param urlPath an URL path or any other string containing the forward slash '/' as directory separator.
     *
     * @return a path string with all occurrences of '/'
     *
     * @throws IllegalArgumentException if the given parameter is <code>null</code>.
     */
    public static String convertToLocalPath(String urlPath) {
        Guardian.assertNotNull("urlPath", urlPath);
        if (File.separatorChar != _URL_DIR_SEPARATOR_CHAR
                && urlPath.indexOf(_URL_DIR_SEPARATOR_CHAR) >= 0) {
            return urlPath.replace(_URL_DIR_SEPARATOR_CHAR,
                                   File.separatorChar);
        }
        return urlPath;
    }

    /**
     * Deletes the directory <code>treeRoot</code> and all the content recursively.
     *
     * @param treeRoot directory to be deleted
     */
    public static void deleteFileTree(File treeRoot) {
        Guardian.assertNotNull("treeRoot", treeRoot);

        File[] files = treeRoot.listFiles();
        if (files != null) {
            for (File file : files) {
                if (file.isDirectory()) {
                    deleteFileTree(file);
                } else {
                    file.delete();
                }
            }
        }
        treeRoot.delete();
    }

    /**
     * Creates a (more) human readable exception message text for the given exception. This method should be used when
     * exception messages are to be presented to the user in a GUI.
     * <p/>
     * <p>Currently the only modifications are<br> 1. the first letter is turned into upper case <br> 2. the message is
     * suffixed with a dot ('.') character.
     *
     * @param e the exception
     *
     * @return a modified message text, or <code>null</code> if <code>e</code> was null.
     */
    public static String createHumanReadableExceptionMessage(final Exception e) {
        if (e == null) {
            return null;
        }
        String message = e.getMessage();
        if (message != null && message.length() > 0) {
            final StringBuffer sb = new StringBuffer();
            sb.append(Character.toUpperCase(message.charAt(0)));
            sb.append(message.substring(1));
            String[] punctuators = new String[]{".", ",", "!", "?", ";", ":"};
            boolean punctuatorFound = false;
            for (String punctuator : punctuators) {
                if (message.endsWith(punctuator)) {
                    punctuatorFound = true;
                    break;
                }
            }
            if (!punctuatorFound) {
                sb.append('.');
            }
            message = sb.toString();
        } else {
            message = "No message text available.";
        }
        return message;
    }

    /**
     * Copies the given text to the system clipboard.
     *
     * @param text the text to copy
     */
    public static void copyToClipboard(final String text) {
        StringSelection selection = new StringSelection(text == null ? "" : text);
        Clipboard clipboard = Toolkit.getDefaultToolkit().getSystemClipboard();
        if (clipboard != null) {
            clipboard.setContents(selection, selection);
        } else {
            BeamLogManager.getSystemLogger().severe("failed to obtain clipboard instance");
        }
    }

    /**
     * Copies the given image to the system clipboard.
     *
     * @param image the image to copy
     */
    public static void copyToClipboard(final Image image) {
        ImageSelection selection = new ImageSelection(image);
        Clipboard clipboard = Toolkit.getDefaultToolkit().getSystemClipboard();
        if (clipboard != null) {
            clipboard.setContents(selection, null);
        } else {
            BeamLogManager.getSystemLogger().severe("failed to obtain clipboard instance");
        }
    }

    public static boolean isRunningOnMacOS() {

        String osName = System.getProperty("os.name");
        if ("Mac OS X".equalsIgnoreCase(osName)) {
            return true;
        }

        final String macOsSpecificPropertyKey = "mrj.version";
        final String systemLafName = UIManager.getSystemLookAndFeelClassName();
        final String currentLafName = UIManager.getLookAndFeel().getClass().getName();

        return System.getProperty(macOsSpecificPropertyKey) != null
                && systemLafName.equals(currentLafName);
    }

    /**
     * Loads services from all <code>META-INF/services/</code> resources.
     *
     * @param serviceType the type of the service to be loaded.
     *
     * @return the services of type <code>serviceType</code> found.
     */
    public static <S> Iterable<S> loadServices(Class<S> serviceType) {
        return ServiceLoader.load(serviceType);
    }

    /**
     * Loads services from all <code>META-INF/services/</code> resources.
     *
     * @param serviceType the type of the service to be loaded.
     * @param classLoader the class loader.
     *
     * @return the services of type <code>serviceType</code> found.
     */
    public static <S> Iterable<S> loadServices(Class<S> serviceType, ClassLoader classLoader) {
        return ServiceLoader.load(serviceType, classLoader);
    }

    public static String getBuildNumber() {
        // todo - in BEAM 3.x org.esa.beam.resources.bundles.build resource has been used. 
        // todo - use application.properties with version ID set by Maven (resource Filter!)
        return "1";
    }

    public static int getLogLevel(String logLevelStr) {
        int logLevel = LL_INFO;
        if (LLS_DEBUG.equalsIgnoreCase(logLevelStr)) {
            logLevel = LL_DEBUG;
        } else if (LLS_INFO.equalsIgnoreCase(logLevelStr)) {
            logLevel = LL_INFO;
        } else if (LLS_ERROR.equalsIgnoreCase(logLevelStr)) {
            logLevel = LL_ERROR;
        } else if (LLS_WARNING.equalsIgnoreCase(logLevelStr)) {
            logLevel = LL_WARNING;
        }
        return logLevel;
    }

    public static Class<?> loadHdf4Lib(Class<?> callerClass) {
        try {
            return Class.forName(_H4_CLASS_NAME, true, callerClass.getClassLoader());
        } catch (Throwable error) {
            BeamLogManager.getSystemLogger().warning(MessageFormat.format("{0}: HDF-4 library not available: {1}: {2}", callerClass, error.getClass(), error.getMessage()));
            return null;
        }
    }

    public static Class<?> loadHdf5Lib(Class<?> callerClass) {
        try {
            return Class.forName(_H5_CLASS_NAME, true, callerClass.getClassLoader());
        } catch (Throwable error) {
            BeamLogManager.getSystemLogger().warning(MessageFormat.format("{0}: HDF-5 library not available: {1}: {2}", callerClass, error.getClass(), error.getMessage()));
            return null;
        }
    }

    /**
     * This class is used to hold an image while on the clipboard.
     */
    public static class ImageSelection implements Transferable {

        private Image _image;

        public ImageSelection(Image image) {
            _image = image;
        }

        // Returns supported flavors
        public DataFlavor[] getTransferDataFlavors() {
            return new DataFlavor[]{DataFlavor.imageFlavor};
        }

        // Returns true if flavor is supported
        public boolean isDataFlavorSupported(DataFlavor flavor) {
            return DataFlavor.imageFlavor.equals(flavor);
        }

        // Returns image
        public Object getTransferData(DataFlavor flavor) throws UnsupportedFlavorException,
                IOException {
            if (!DataFlavor.imageFlavor.equals(flavor)) {
                throw new UnsupportedFlavorException(flavor);
            }
            return _image;
        }
    }

}
