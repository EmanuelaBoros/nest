/*
 * $Id: MaskViewerToolView.java,v 1.2 2009-12-04 19:06:45 lveci Exp $
 *
 * Copyright (C) 2009 by Brockmann Consult (info@brockmann-consult.de)
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
package org.esa.beam.visat.toolviews.mask;

import org.esa.beam.framework.ui.application.support.AbstractToolView;

public class MaskViewerToolView extends MaskToolView {

    public static final String ID = MaskViewerToolView.class.getName();

    @Override
    protected MaskForm createMaskForm(AbstractToolView maskToolView) {
        return new MaskViewerForm();
    }
}