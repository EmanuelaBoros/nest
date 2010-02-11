/*
 * $Id: InternalFrameLayoutManager.java,v 1.2 2010-02-11 18:32:53 lveci Exp $
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
package com.bc.swing.desktop;

import javax.swing.JDesktopPane;
import javax.swing.JInternalFrame;

/**
 * @author Norman Fomferra (norman.fomferra@brockmann-consult.de)
 * @version $Revision: 1.2 $ $Date: 2010-02-11 18:32:53 $
 */
public interface InternalFrameLayoutManager {
    void moveFrameToVisible(JDesktopPane desktopPane, JInternalFrame frame);
    void cascadeFrames(JDesktopPane desktopPane, JInternalFrame[] frames);
    void tileFramesEvenly(JDesktopPane desktopPane, JInternalFrame[] frames);
    void tileFramesHorizontally(JDesktopPane desktopPane, JInternalFrame[] frames);
    void tileFramesVertically(JDesktopPane desktopPane, JInternalFrame[] frames);
}
