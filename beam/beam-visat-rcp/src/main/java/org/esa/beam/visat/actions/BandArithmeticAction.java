package org.esa.beam.visat.actions;

import org.esa.beam.framework.datamodel.Product;
import org.esa.beam.framework.datamodel.ProductNodeList;
import org.esa.beam.framework.ui.command.CommandEvent;
import org.esa.beam.framework.ui.command.ExecCommand;
import org.esa.beam.visat.VisatApp;
import org.esa.beam.visat.dialogs.BandArithmetikDialog;


/**
 * VISAT's band arithmetic feature.
 *
 * @author Norman Fomferra
 * @version $Revision: 1.2 $ $Date: 2009-04-29 19:35:16 $
 */
public class BandArithmeticAction extends ExecCommand {

    @Override
    public void actionPerformed(final CommandEvent event) {
        openBandArithmeticDialog(VisatApp.getApp(), event.getCommand().getHelpId());
    }

    @Override
    public void updateState(final CommandEvent event) {
        final int n = VisatApp.getApp().getProductManager().getProductCount();
        setEnabled(n > 0);
    }

   
    private void openBandArithmeticDialog(final VisatApp visatApp, final String helpId) {

        final Product[] prods = visatApp.getProductManager().getProducts();
        final ProductNodeList<Product> products = new ProductNodeList<Product>();
        for (Product prod : prods) {
            products.add(prod);
        }
        BandArithmetikDialog bandArithmetikDialog = new BandArithmetikDialog(visatApp,
                                                         visatApp.getSelectedProduct(),
                                                         products,
                                                         helpId);
        bandArithmetikDialog.show();
    }
}
