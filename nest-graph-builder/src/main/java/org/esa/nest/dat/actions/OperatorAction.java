package org.esa.nest.dat.actions;

import org.esa.beam.visat.actions.DefaultOperatorAction;
import org.esa.beam.framework.ui.ModelessDialog;
import org.esa.nest.dat.dialogs.NestSingleTargetProductDialog;
import org.esa.nest.util.ResourceUtils;
import com.bc.ceres.core.runtime.ConfigurationElement;
import com.bc.ceres.core.CoreException;

import javax.swing.*;

/**
 * <p>An action which creates a default dialog for an operator given by the
 * action property action property {@code operatorName}.</p>
 * <p>Optionally the dialog title can be set via the {@code dialogTitle} property and
 * the ID of the help page can be given using the {@code helpId} property. If not given the
 * name of the operator will be used instead. Also optional the
 * file name suffix for the target product can be given via the {@code targetProductNameSuffix} property.</p>
 */
public class OperatorAction extends DefaultOperatorAction {
    private String iconName;

    @Override
    public void configure(ConfigurationElement config) throws CoreException {
        super.configure(config);
        iconName = getConfigString(config, "icon");
    }

    @Override
    protected ModelessDialog createOperatorDialog() {
        final NestSingleTargetProductDialog productDialog = new NestSingleTargetProductDialog(operatorName,
                                                    getAppContext(),  dialogTitle, getHelpId());
        if (targetProductNameSuffix != null) {
            productDialog.setTargetProductNameSuffix(targetProductNameSuffix);
        }
        addIcon(productDialog);
        return productDialog;
    }

    protected void addIcon(final NestSingleTargetProductDialog productDialog) {
        if(iconName == null) {
            productDialog.setIcon(ResourceUtils.nestIcon);
        } else if(iconName.equals("rstbIcon")) {
            productDialog.setIcon(ResourceUtils.rstbIcon); 
        } else {
            final ImageIcon icon = ResourceUtils.LoadIcon(iconName);
            if(icon != null)
                productDialog.setIcon(icon);
        }
    }
}
