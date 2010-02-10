package org.esa.beam.visat.toolviews.placemark;

import com.bc.ceres.core.Assert;
import org.esa.beam.framework.datamodel.Placemark;
import org.esa.beam.framework.datamodel.PlacemarkDescriptor;
import org.esa.beam.framework.datamodel.Product;
import org.esa.beam.framework.datamodel.ProductNodeGroup;

/**
 * Created by Marco Peters.
 *
 * @author Marco Peters
 * @version $Revision: 1.2 $ $Date: 2010-02-10 16:20:37 $
 */
public class PlacemarkNameFactory {

    private PlacemarkNameFactory() {
    }

    public static String[] createUniqueNameAndLabel(PlacemarkDescriptor placemarkDescriptor, Product product) {
        ProductNodeGroup<Placemark> placemarkGroup = placemarkDescriptor.getPlacemarkGroup(product);
        return createUniqueNameAndLabel(placemarkDescriptor, placemarkGroup);
    }

    public static String[] createUniqueNameAndLabel(PlacemarkDescriptor placemarkDescriptor, ProductNodeGroup<Placemark> group) {
        int pinNumber = group.getNodeCount() + 1;
        String name = PlacemarkNameFactory.createName(placemarkDescriptor, pinNumber);
        while (group.get(name) != null) {
            name = PlacemarkNameFactory.createName(placemarkDescriptor, ++pinNumber);
        }
        final String label = PlacemarkNameFactory.createLabel(placemarkDescriptor, pinNumber, true);
        return new String[]{name, label};
    }

    public static String createLabel(PlacemarkDescriptor placemarkDescriptor, int pinNumber, boolean firstLetterIsUpperCase) {
        Assert.argument(placemarkDescriptor.getRoleLabel().length() > 0, "placemarkDescriptor.getRoleLabel()");
        String name = placemarkDescriptor.getRoleLabel();
        if (firstLetterIsUpperCase) {
            name = name.substring(0, 1).toUpperCase() + name.substring(1);
        }
        return name + " " + pinNumber;
    }

    public static String createName(PlacemarkDescriptor placemarkDescriptor, int pinNumber) {
        return placemarkDescriptor.getRoleName() + "_" + pinNumber;
    }

}
