package com.bc.ceres.swing.binding.internal;

import com.bc.ceres.binding.PropertyDescriptor;
import com.bc.ceres.swing.binding.ComponentAdapter;

import javax.swing.JComponent;
import javax.swing.JSpinner;
import javax.swing.SpinnerListModel;
import javax.swing.SpinnerNumberModel;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

/**
 * A binding for a {@link JSpinner} component.
 *
 * @author Norman Fomferra
 * @version $Revision: 1.1 $ $Date: 2010-02-10 19:57:11 $
 * @since BEAM 4.2
 */
public class SpinnerAdapter extends ComponentAdapter implements ChangeListener {

    final JSpinner spinner;

    public SpinnerAdapter(JSpinner spinner) {
        this.spinner = spinner;
    }

    @Override
    public JComponent[] getComponents() {
        return new JComponent[]{spinner};
    }

    @Override
    public void bindComponents() {
        updateSpinnerModel();
        spinner.addChangeListener(this);
    }

    @Override
    public void unbindComponents() {
        spinner.removeChangeListener(this);
    }

    @Override
    public void adjustComponents() {
        Object value = getBinding().getPropertyValue();
        spinner.setValue(value);
    }

    @Override
    public void stateChanged(ChangeEvent evt) {
        getBinding().setPropertyValue(spinner.getValue());
    }

    private void updateSpinnerModel() {
        PropertyDescriptor propertyDescriptor = getBinding().getContext().getPropertySet().getDescriptor(getBinding().getPropertyName());
        if (propertyDescriptor.getValueRange() != null) {
            Class<?> type = propertyDescriptor.getType();

            if (Number.class.isAssignableFrom(type)) {
                Number defaultValue = (Number) propertyDescriptor.getDefaultValue(); // todo - why not the current value? (mp,nf - 18.02.2008)
                double min = propertyDescriptor.getValueRange().getMin();
                double max = propertyDescriptor.getValueRange().getMax();
                // todo - get step size from interval

                if (type == Byte.class) {
                    spinner.setModel(new SpinnerNumberModel(defaultValue, (byte) min, (byte) max, 1));
                } else if (type == Short.class) {
                    spinner.setModel(new SpinnerNumberModel(defaultValue, (short) min, (short) max, 1));
                } else if (type == Integer.class) {
                    spinner.setModel(new SpinnerNumberModel(defaultValue, (int) min, (int) max, 1));
                } else if (type == Long.class) {
                    spinner.setModel(new SpinnerNumberModel(defaultValue, (long) min, (long) max, 1));
                } else if (type == Float.class) {
                    spinner.setModel(new SpinnerNumberModel(defaultValue, (float) min, (float) max, 1));
                } else {
                    spinner.setModel(new SpinnerNumberModel(defaultValue, min, max, 1));
                }
            }
        } else if (propertyDescriptor.getValueSet() != null) {
            spinner.setModel(new SpinnerListModel(propertyDescriptor.getValueSet().getItems()));
        }
    }
}
