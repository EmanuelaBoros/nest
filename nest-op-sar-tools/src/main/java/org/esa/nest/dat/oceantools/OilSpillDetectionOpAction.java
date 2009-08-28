package org.esa.nest.dat.oceantools;

import org.esa.beam.framework.gpf.ui.DefaultSingleTargetProductDialog;
import org.esa.beam.framework.ui.ModelessDialog;
import org.esa.beam.framework.ui.command.CommandEvent;
import org.esa.beam.visat.actions.AbstractVisatAction;

/**
 * Oil Spill Detection action.
 *
 */
public class OilSpillDetectionOpAction extends AbstractVisatAction {

    private ModelessDialog dialog;

    @Override
    public void actionPerformed(CommandEvent event) {

        if (dialog == null) {
            dialog = new DefaultSingleTargetProductDialog("Oil-Spill-Detection", getAppContext(), "Oil-Spill-Detection", getHelpId());
        }
        dialog.show();

    }

    @Override
    public void updateState(final CommandEvent event) {
        setEnabled(true);
    }
}