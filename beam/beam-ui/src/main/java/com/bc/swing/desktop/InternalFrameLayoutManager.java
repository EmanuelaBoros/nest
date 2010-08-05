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
package com.bc.swing.desktop;

import javax.swing.JDesktopPane;
import javax.swing.JInternalFrame;

/**
 * @author Norman Fomferra (norman.fomferra@brockmann-consult.de)
 * @version $Revision: 1.3 $ $Date: 2010-08-05 17:00:54 $
 */
public interface InternalFrameLayoutManager {
    void moveFrameToVisible(JDesktopPane desktopPane, JInternalFrame frame);
    void cascadeFrames(JDesktopPane desktopPane, JInternalFrame[] frames);
    void tileFramesEvenly(JDesktopPane desktopPane, JInternalFrame[] frames);
    void tileFramesHorizontally(JDesktopPane desktopPane, JInternalFrame[] frames);
    void tileFramesVertically(JDesktopPane desktopPane, JInternalFrame[] frames);
}
