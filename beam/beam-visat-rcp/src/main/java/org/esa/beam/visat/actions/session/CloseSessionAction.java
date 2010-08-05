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
package org.esa.beam.visat.actions.session;

import org.esa.beam.framework.ui.command.CommandEvent;
import org.esa.beam.framework.ui.command.ExecCommand;
import org.esa.beam.visat.VisatApp;


/**
 * Closes a VISAT session.
 *
 * @author Norman Fomferra
 * @version $Revision: 1.7 $ $Date: 2010-08-05 17:00:55 $
 * @since BEAM 4.6
 */
public class CloseSessionAction extends ExecCommand {

    public static final String ID = "closeSession";

    @Override
    public void actionPerformed(final CommandEvent event) {
        final VisatApp app = VisatApp.getApp();
        app.closeAllProducts();
        app.setSessionFile(null);
    }

    @Override
    public void updateState(final CommandEvent event) {
        final VisatApp app = VisatApp.getApp();
        setEnabled(app.getProductManager().getProductCount() > 0 || app.getSessionFile() != null);
    }
}