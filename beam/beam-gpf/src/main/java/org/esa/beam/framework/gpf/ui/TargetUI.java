package org.esa.beam.framework.gpf.ui;

import org.esa.beam.framework.ui.AppContext;
import org.esa.beam.framework.ui.BasicApp;
import org.esa.beam.util.SystemUtils;

import javax.swing.*;
import java.util.Map;
import java.io.File;

/**
 * Writer OperatorUI
 */
public class TargetUI extends BaseOperatorUI {

    TargetProductSelector targetProductSelector = null;
    private static final String FILE_PARAMETER = "file";
    private AppContext appContext;

    @Override
    public JComponent CreateOpTab(String operatorName, Map<String, Object> parameterMap, AppContext appContext) {

        paramMap = parameterMap;
        targetProductSelector = new TargetProductSelector();
        this.appContext = appContext;

        File saveDir = null;
        final Object value = paramMap.get(FILE_PARAMETER);
        if(value != null) {
            final File file = (File)value;
            saveDir = file.getParentFile();
        }

        if(saveDir == null) {
            final String homeDirPath = SystemUtils.getUserHomeDir().getPath();
            final String savePath = appContext.getPreferences().getPropertyString(BasicApp.PROPERTY_KEY_APP_LAST_SAVE_DIR, homeDirPath);
            saveDir = new File(savePath);
        }
        targetProductSelector.getModel().setProductDir(saveDir);
        targetProductSelector.getOpenInAppCheckBox().setText("Open in " + appContext.getApplicationName());

        initParameters();
        
        return targetProductSelector.createDefaultPanel();
    }

    @Override
    public void initParameters() {
        assert(paramMap != null);
        final Object value = paramMap.get(FILE_PARAMETER);
        if(value != null) {
            final File file = (File)value;
            targetProductSelector.getProductNameTextField().setText(file.getName());
            targetProductSelector.getModel().setProductName(file.getName());
        }
    }

    @Override
    public UIValidation validateParameters() {

        final String productName = targetProductSelector.getModel().getProductName();
        if(productName == null || productName.isEmpty())
            return new UIValidation(UIValidation.State.ERROR, "Target file not specified");
        final File file = targetProductSelector.getModel().getProductFile();
        if(file == null)
            return new UIValidation(UIValidation.State.ERROR, "Target file not specified");

        final String productDir = targetProductSelector.getModel().getProductDir().getAbsolutePath();
        appContext.getPreferences().setPropertyString(BasicApp.PROPERTY_KEY_APP_LAST_SAVE_DIR, productDir);

        return new UIValidation(UIValidation.State.OK, "");
    }

    @Override
    public void updateParameters() {

        if(targetProductSelector.getModel().getProductName() != null) {
            paramMap.put("file", targetProductSelector.getModel().getProductFile());
            paramMap.put("formatName", targetProductSelector.getModel().getFormatName());
        }
    }
}
