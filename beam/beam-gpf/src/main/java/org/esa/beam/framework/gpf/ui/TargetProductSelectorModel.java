package org.esa.beam.framework.gpf.ui;

import com.bc.ceres.binding.PropertyContainer;
import com.bc.ceres.binding.PropertyDescriptor;
import com.bc.ceres.binding.validators.NotEmptyValidator;
import com.bc.ceres.binding.validators.NotNullValidator;
import org.esa.beam.framework.dataio.ProductIO;
import org.esa.beam.framework.dataio.ProductIOPlugInManager;
import org.esa.beam.framework.dataio.ProductWriterPlugIn;
import org.esa.beam.util.StringUtils;

import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.File;
import java.util.Iterator;

/**
 * WARNING: This class belongs to a preliminary API and may change in future releases.
 *
 * @author Ralf Quast
 * @version $Revision: 1.2 $ $Date: 2009-11-04 17:04:32 $
 */
@SuppressWarnings({"UnusedDeclaration"})
public class TargetProductSelectorModel {

    // used by object binding
    private String productName;
    private boolean saveToFileSelected;
    private boolean openInAppSelected;
    private File productDir;
    private String formatName;
    private String[] formatNames = { ProductIO.DEFAULT_FORMAT_NAME };

    private final PropertyContainer propertyContainer;

    TargetProductSelectorModel() {
        propertyContainer = PropertyContainer.createObjectBacked(this);
        propertyContainer.addPropertyChangeListener("saveToFileSelected", new PropertyChangeListener() {
            public void propertyChange(PropertyChangeEvent evt) {
                if (!(Boolean) evt.getNewValue()) {
                    setOpenInAppSelected(true);
                }
            }
        });
        propertyContainer.addPropertyChangeListener("openInAppSelected", new PropertyChangeListener() {
            public void propertyChange(PropertyChangeEvent evt) {
                if (!(Boolean) evt.getNewValue()) {
                    setSaveToFileSelected(true);
                }
            }
        });
        PropertyDescriptor productNameDescriptor = propertyContainer.getDescriptor("productName");
        productNameDescriptor.setValidator(new NotEmptyValidator());
        productNameDescriptor.setDisplayName("target product name");

        PropertyDescriptor productDirDescriptor = propertyContainer.getDescriptor("productDir");
        productDirDescriptor.setValidator(new NotNullValidator());
        productDirDescriptor.setDisplayName("target product directory");

        setOpenInAppSelected(true);
        setSaveToFileSelected(true);
        //formatNames = ProductIOPlugInManager.getInstance().getAllProductWriterFormatStrings();
        if (StringUtils.contains(formatNames, ProductIO.DEFAULT_FORMAT_NAME)) {
            setFormatName(ProductIO.DEFAULT_FORMAT_NAME);
        } else {
            setFormatName(formatNames[0]);
        }
    }

    public String getProductName() {
        return productName;
    }

    public boolean isSaveToFileSelected() {
        return saveToFileSelected;
    }

    public boolean isOpenInAppSelected() {
        return openInAppSelected;
    }

    public File getProductDir() {
        return productDir;
    }

    public File getProductFile() {
        return new File(productDir, getProductFileName());
    }

    String getProductFileName() {
        String productFileName = productName;
        Iterator<ProductWriterPlugIn> iterator = ProductIOPlugInManager.getInstance().getWriterPlugIns(formatName);
        if (iterator.hasNext()) {
            final ProductWriterPlugIn writerPlugIn = iterator.next();

            boolean ok = false;
            for (String extension : writerPlugIn.getDefaultFileExtensions()) {
                if (productFileName.endsWith(extension)) {
                    ok = true;
                    break;
                }
            }
            if (!ok) {
                productFileName = productFileName.concat(writerPlugIn.getDefaultFileExtensions()[0]);
            }
        }
        return productFileName;
    }

    public String getFormatName() {
        return formatName;
    }

    public String[] getFormatNames() {
        return formatNames;
    }

    public void setProductName(String productName) {
        setValueContainerValue("productName", productName);
    }

    public void setSaveToFileSelected(boolean saveToFileSelected) {
        setValueContainerValue("saveToFileSelected", saveToFileSelected);
    }

    public void setOpenInAppSelected(boolean openInAppSelected) {
        setValueContainerValue("openInAppSelected", openInAppSelected);
    }

    public void setProductDir(File productDir) {
        setValueContainerValue("productDir", productDir);
    }

    public void setFormatName(String formatName) {
        setValueContainerValue("formatName", formatName);
    }

    public PropertyContainer getValueContainer() {
        return propertyContainer;
    }

    private void setValueContainerValue(String name, Object value) {
        propertyContainer.setValue(name, value);
    }

}
