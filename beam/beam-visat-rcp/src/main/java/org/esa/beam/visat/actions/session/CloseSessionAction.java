/*
 * $Id: CloseSessionAction.java,v 1.3 2009-05-12 12:56:42 lveci Exp $
 *
 * Copyright (C) 2002 by Brockmann Consult (info@brockmann-consult.de)
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation. This program is distributed in the hope it will
 * be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */
package org.esa.beam.visat.actions.session;

import org.esa.beam.framework.ui.command.CommandEvent;
import org.esa.beam.framework.ui.command.ExecCommand;
import org.esa.beam.visat.VisatApp;


/**
 * Closes a VISAT session.
 *
 * @author Norman Fomferra
 * @version $Revision: 1.3 $ $Date: 2009-05-12 12:56:42 $
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
        setEnabled(VisatApp.getApp().getProductManager().getProductCount() > 0);
    }
}