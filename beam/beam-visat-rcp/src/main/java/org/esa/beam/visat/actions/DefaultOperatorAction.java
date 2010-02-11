package org.esa.beam.visat.actions;

import com.bc.ceres.core.CoreException;
import com.bc.ceres.core.runtime.ConfigurationElement;

import org.esa.beam.framework.gpf.ui.DefaultSingleTargetProductDialog;
import org.esa.beam.framework.ui.ModelessDialog;
import org.esa.beam.framework.ui.command.CommandEvent;


/**
 * <p><b>WARNING:</b> This class belongs to a preliminary API and may change in future releases.<p/>
 * 
 * <p>An action which creates a default dialog for an operator given by the
 * action property action property {@code operatorName}.</p>
 * <p>Optionally the dialog title can be set via the {@code dialogTitle} property and
 * the ID of the help page can be given using the {@code helpId} property. If not given the
 * name of the operator will be used instead. Also optional the 
 * file name suffix for the target product can be given via the {@code nameSuffix} property.</p>
 *
 * @author Norman Fomferra
 * @author Marco Zuehlke
 * @version $Revision: 1.1 $ $Date: 2010-02-11 17:02:24 $
 */
public class DefaultOperatorAction extends AbstractVisatAction {

    private ModelessDialog dialog;
    private String operatorName;
    private String dialogTitle;
    private String nameSuffix;

    @Override
    public void actionPerformed(CommandEvent event) {
      if (dialog == null) {
            dialog = createOperatorDialog();
        }
        dialog.show();
    }
    
    @Override
    public void configure(ConfigurationElement config) throws CoreException {
        operatorName = getConfigString(config, "operatorName");
        if (operatorName == null) {
            throw new CoreException("Missing DefaultOperatorAction property 'operatorName'.");
        }
        dialogTitle = getValue(config, "dialogTitle", operatorName);
        nameSuffix = getConfigString(config, "nameSuffix");
        super.configure(config);
    }

    protected ModelessDialog createOperatorDialog() {
        DefaultSingleTargetProductDialog productDialog = new DefaultSingleTargetProductDialog(operatorName, getAppContext(),
                                                    dialogTitle, getHelpId());
        if (nameSuffix != null) {
            productDialog.setTargetProductNameSuffix(nameSuffix);
        }
        return productDialog;
    }
}
