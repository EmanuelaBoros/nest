/*
 * Copyright (C) 2010 Brockmann Consult GmbH (info@brockmann-consult.de)
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation; either version 3 of the License, or (at your option)
 * any later version.
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, see http://www.gnu.org/licenses/
 */

package org.esa.beam.framework.gpf.main;

import org.esa.beam.framework.dataio.ProductIO;
import org.esa.beam.framework.dataio.ProductIOPlugInManager;
import org.esa.beam.framework.dataio.ProductWriterPlugIn;
import org.esa.beam.framework.gpf.GPF;
import org.esa.beam.util.io.FileUtils;

import java.text.MessageFormat;
import java.util.Iterator;
import java.util.SortedMap;
import java.util.TreeMap;

/**
 * The common command-line for GPF.
 */
class CommandLineArgs {

    private static final int K = 1024;
    private static final int M = K * 1024;
    private static final int G = M * 1024;

    private String[] args;
    private String operatorName;
    private String graphFilepath;
    private String targetFilepath;
    private TreeMap<String, String> parameterMap;
    private TreeMap<String, String> sourceFilepathMap;
    private String targetFormatName;
    private String parameterFilepath;
    private String fileListPath = null;
    private TreeMap<String, String> targetFilepathMap;
    private boolean helpRequested;
    private boolean printAllHelp = false;
    private boolean stackTraceDump;
    private boolean clearCacheAfterRowWrite;
    private long tileCacheCapacity;

    private int tileSchedulerParallelism;

    public CommandLineArgs(String[] args) {
        this.args = args.clone();
        if (this.args.length == 0) {
            helpRequested = true;
        }

        sourceFilepathMap = new TreeMap<String, String>();
        targetFilepathMap = new TreeMap<String, String>();
        parameterMap = new TreeMap<String, String>();
        tileCacheCapacity = CommandLineTool.DEFAULT_TILE_CACHE_SIZE_IN_M * M;
        tileSchedulerParallelism = CommandLineTool.DEFAULT_TILE_SCHEDULER_PARALLELSIM;

        // look for "-e" early do enable verbose error reports 
        for (String arg : this.args) {
            if (arg.equals("-e")) {
                stackTraceDump = true;
            }
        }
    }

    public void parseArguments() throws Exception {
        int argCount = 0;
        for (int i = 0; i < this.args.length; i++) {
            String arg = this.args[i];
            if (arg.startsWith("-")) {
                if (arg.startsWith("-P")) {
                    String[] pair = parseNameValuePair(arg);
                    parameterMap.put(pair[0], pair[1]);
                } else if (arg.startsWith("-S")) {
                    String[] pair = parseNameValuePair(arg);
                    sourceFilepathMap.put(pair[0], pair[1]);
                } else if (arg.startsWith("-T")) {
                    String[] pair = parseNameValuePair(arg);
                    targetFilepathMap.put(pair[0], pair[1]);
                } else if (arg.equals("-h")) {
                    helpRequested = true;
                } else if (arg.equals("-printHelp")) {
                    printAllHelp = true;
                } else if (arg.equals("-x")) {
                    clearCacheAfterRowWrite = true;
                } else if (arg.equals("-e")) {
                    // already parsed
                } else if (arg.equals("-t")) {
                    targetFilepath = parseOptionArgument(arg, i);
                    i++;
                } else if (arg.equals("-f")) {
                    targetFormatName = parseOptionArgument(arg, i);
                    i++;
                } else if (arg.equals("-p")) {
                    parameterFilepath = parseOptionArgument(arg, i);
                    i++;
                } else if (arg.equals("-q")) {
                    tileSchedulerParallelism = parseOptionArgumentInt(arg, i);
                    i++;
                } else if (arg.equals("-c")) {
                    tileCacheCapacity = parseOptionArgumentBytes(arg, i);
                    i++;
                } else if (arg.equals("-fileListPath")) {
                    fileListPath = parseOptionArgument(arg, i);
                    i++;
                } else {
                    throw error("Unknown option '" + arg + "'");
                }
            } else {
                if (argCount == 0) {
                    if (arg.endsWith(".xml") || arg.endsWith(".XML") || arg.contains("/") || arg.contains("\\")) {
                        graphFilepath = arg;
                    } else {
                        operatorName = arg;
                    }
                } else {
                    int index = argCount - 1;
                    if (index == 0) {
                        sourceFilepathMap.put(GPF.SOURCE_PRODUCT_FIELD_NAME, arg);
                    }
                    sourceFilepathMap.put(GPF.SOURCE_PRODUCT_FIELD_NAME + "." + (index + 1), arg);
                    // kept for backward compatibility
                    // since BEAM 4.9 the pattern above is preferred
                    sourceFilepathMap.put(GPF.SOURCE_PRODUCT_FIELD_NAME + (index + 1), arg);
                }
                argCount++;
            }
        }

        if (operatorName == null && graphFilepath == null && !helpRequested && !printAllHelp) {
            throw error("Either operator name or graph XML file must be given");
        }
        if (graphFilepath == null && !targetFilepathMap.isEmpty()) {
            throw error("Defined target products only valid for graph XML");
        }
        if (targetFilepath == null && targetFilepathMap.isEmpty()) {
            targetFilepath = CommandLineTool.DEFAULT_TARGET_FILEPATH;
        }
        if (targetFormatName == null && targetFilepath != null) {
            final String extension = FileUtils.getExtension(targetFilepath);
            if (extension == null || extension.isEmpty()) {
                targetFormatName = ProductIO.DEFAULT_FORMAT_NAME;
                // todo - decide if to append extension or not  (nf - 29.10.2007)
            } else {
                targetFormatName = detectWriterFormat(extension);
                if (targetFormatName == null) {
                    throw error("Output format unknown");
                }
            }
        }
    }

    public String[] getArgs() {
        return args;
    }

    public String getOperatorName() {
        return operatorName;
    }

    public String getGraphFilepath() {
        return graphFilepath;
    }

    public String getTargetFilepath() {
        return targetFilepath;
    }

    public String getTargetFormatName() {
        return targetFormatName;
    }

    public String getParameterFilepath() {
        return parameterFilepath;
    }

    public String getFileListPath() {
        return fileListPath;
    }

    public long getTileCacheCapacity() {
        return tileCacheCapacity;
    }

    public int getTileSchedulerParallelism() {
        return tileSchedulerParallelism;
    }

    public boolean isClearCacheAfterRowWrite() {
        return clearCacheAfterRowWrite;
    }

    public SortedMap<String, String> getParameterMap() {
        return parameterMap;
    }

    public SortedMap<String, String> getSourceFilepathMap() {
        return sourceFilepathMap;
    }

    public SortedMap<String, String> getTargetFilepathMap() {
        return targetFilepathMap;
    }

    public boolean isHelpRequested() {
        return helpRequested;
    }

    public boolean isPrintAllHelpRequested() {
        return printAllHelp;
    }

    public boolean isStackTraceDump() {
        return stackTraceDump;
    }

    private String parseOptionArgument(String arg, int index) throws Exception {
        if (index < args.length - 1) {
            return args[index + 1];
        } else {
            throw error("Missing argument for option '" + arg + "'");
        }
    }

    private int parseOptionArgumentInt(String arg, int index) throws Exception {
       String valueString = parseOptionArgument(arg, index);
       return Integer.parseInt(valueString);
    }

    private long parseOptionArgumentBytes(String arg, int index) throws Exception {
        String valueString = parseOptionArgument(arg, index);
        long factor = 1;
        if (valueString.toUpperCase().endsWith("K")) {
            factor = K;
            valueString = valueString.substring(0, valueString.length() - 1);
        } else if (valueString.toUpperCase().endsWith("M")) {
            factor = M;
            valueString = valueString.substring(0, valueString.length() - 1);
        } else if (valueString.toUpperCase().endsWith("G")) {
            factor = G;
            valueString = valueString.substring(0, valueString.length() - 1);
        }

        long value = Long.parseLong(valueString);
        if (value < 0L) {
            throw error(MessageFormat.format("Value for ''{0}'' must not be negative", arg));
        }
        return factor * value;
    }

    private String[] parseNameValuePair(String arg) throws Exception {
        int pos = arg.indexOf('=');
        if (pos == -1) {
            throw error("Missing '=' in '" + arg + "'");
        }
        String name = arg.substring(2, pos).trim();
        if (name.isEmpty()) {
            throw error("Empty identifier in '" + arg + "'");
        }
        String value = arg.substring(pos + 1).trim();
        return new String[]{name, value};
    }

    private String detectWriterFormat(String extension) {
        ProductIOPlugInManager registry = ProductIOPlugInManager.getInstance();
        Iterator<ProductWriterPlugIn> ins = registry.getAllWriterPlugIns();
        while (ins.hasNext()) {
            ProductWriterPlugIn productWriterPlugIn = ins.next();
            String[] strings = productWriterPlugIn.getDefaultFileExtensions();
            for (String string : strings) {
                if (string.equalsIgnoreCase(extension)) {
                    return productWriterPlugIn.getFormatNames()[0];
                }
            }
        }
        return null;
    }

    private static Exception error(String m) {
        return new Exception(m);
    }
}
