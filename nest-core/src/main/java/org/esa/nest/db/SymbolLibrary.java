package org.esa.nest.db;

import org.esa.nest.util.ResourceUtils;
import org.esa.beam.visat.VisatApp;

import javax.swing.*;
import java.io.File;
import java.util.Map;
import java.util.HashMap;
import java.awt.image.BufferedImage;

/**
 * Caches symbols
 */
public class SymbolLibrary {

    private final static File symbolFolder = new File(ResourceUtils.getResFolder(), "symbols");
    private final static String symbolFolderStr = symbolFolder.getAbsolutePath();
    private final static String symbolFolderPlaceHolder = "{$SYMBOL_FOLDER}";
    private final static File defaultSymbolFile = new File(symbolFolder, "flag_blue.png");

    private static SymbolLibrary _instance = null;
    private static final Map<File, SymbolNode> fileMap = new HashMap<File, SymbolNode>(20);
    private final SymbolNode rootNode;

    public static SymbolLibrary instance() {
        if(_instance == null) {
            _instance = new SymbolLibrary();
        }
        return _instance;
    }

    private SymbolLibrary() {
        rootNode = new SymbolNode(symbolFolder);
        rootNode.scanSymbolFolder(symbolFolder, fileMap);
    }

    public SymbolNode getRootNode() {
        return rootNode;
    }

    public static SymbolNode getSymbol(final File file) {
        return fileMap.get(file);
    }

    public static ImageIcon getIcon(final File file) {
        return fileMap.get(file).getIcon();
    }

    public static BufferedImage getImage(final File file) {
        return fileMap.get(file).getImage();
    }

    public static SymbolNode getDefaultSymbol() {
        SymbolNode defaultSymbol = fileMap.get(getLastSymbolUsed());
        if(defaultSymbol == null) {
            defaultSymbol = fileMap.get(fileMap.keySet().iterator().next());
        }
        return defaultSymbol;
    }

    private static File getLastSymbolUsed() {
        if(VisatApp.getApp() != null) {
            final String path = VisatApp.getApp().getPreferences().getPropertyString("last_symbol_used", "");
            if(path != null && !path.isEmpty()) {
                return new File(path);
            }
        }
        return defaultSymbolFile;
    }

    public static void saveLastSymbolUsed(final File file) {
        if(VisatApp.getApp() != null && file != null) {
            VisatApp.getApp().getPreferences().setPropertyString("last_symbol_used", file.getAbsolutePath());
        }
    }

    public static String translatePath(final File file) {
        if(file == null) return "";
        String filePath = file.getAbsolutePath();
        if(filePath.startsWith(symbolFolderStr)) {
            filePath = filePath.replace(symbolFolderStr, symbolFolderPlaceHolder);
        }
        return filePath;
    }

    public static File translatePath(String path) {
        if(path == null || path.isEmpty()) return null;
        if(path.startsWith(symbolFolderPlaceHolder)) {
            path = path.replace(symbolFolderPlaceHolder, symbolFolderStr);
        }
        return new File(path);
    }

    public static String clipPath(final File file) {
        String path = file.getAbsolutePath();
        if(path.equals(symbolFolderStr))
            return File.separator;
        else if(path.startsWith(symbolFolderStr)) {
            path = path.replace(symbolFolderStr, "");
        }
        return path;
    }
}
