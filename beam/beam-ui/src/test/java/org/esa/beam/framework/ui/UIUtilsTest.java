/*
 * $Id: UIUtilsTest.java,v 1.1 2009-04-28 14:17:18 lveci Exp $
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

package org.esa.beam.framework.ui;

import junit.framework.TestCase;
import org.esa.beam.framework.param.ParamProperties;
import org.esa.beam.framework.param.Parameter;

import javax.swing.JInternalFrame;
import javax.swing.JSpinner;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.HeadlessException;
import java.awt.Panel;
import java.awt.Rectangle;

public class UIUtilsTest extends TestCase {

    public void testCenterComponent() {
        try {
            Dimension screenSize = java.awt.Toolkit.getDefaultToolkit().getScreenSize();
            Component comp = new Panel();
            Component alignComp = new Panel();

            comp.setBounds(0, 0, 100, 100);
            UIUtils.centerComponent(comp);
            assertEquals(comp.getBounds(), new Rectangle(screenSize.width / 2 - 50, screenSize.height / 2 - 50, 100, 100));

            comp.setBounds(0, 0, 100, 100);
            alignComp.setBounds(100, 100, 200, 200);
            UIUtils.centerComponent(comp, alignComp);
            assertEquals(comp.getBounds(), new Rectangle(150, 150, 100, 100));
        } catch (HeadlessException e) {
            warnHeadless();
        }
    }

    public void testGetScreenSize() {
        try {
            Dimension screenSize = java.awt.Toolkit.getDefaultToolkit().getScreenSize();
            assertNotNull(UIUtils.getScreenSize());
            assertEquals(UIUtils.getScreenSize(), screenSize);
        } catch (HeadlessException e) {
            warnHeadless();
        }
    }

    public void testGetScreenWidth() {
        try {
            Dimension screenSize = java.awt.Toolkit.getDefaultToolkit().getScreenSize();
            assertEquals(UIUtils.getScreenWidth(), screenSize.width);
        } catch (HeadlessException e) {
            warnHeadless();
        }
    }

    public void testGetScreenHeight() {
        try {
            Dimension screenSize = java.awt.Toolkit.getDefaultToolkit().getScreenSize();
            assertEquals(UIUtils.getScreenHeight(), screenSize.height);
        } catch (HeadlessException e) {
            warnHeadless();
        }
    }


    public void testGetUniqueFrameTitle() {
        String title;

        title = UIUtils.getUniqueFrameTitle(
                new JInternalFrame[]{
                        new JInternalFrame("Image1"),
                        new JInternalFrame("Image2"),
                        new JInternalFrame("Image3"),
                }, "Image");
        assertEquals("Image", title);

        title = UIUtils.getUniqueFrameTitle(
                new JInternalFrame[]{
                        new JInternalFrame("Image"),
                        new JInternalFrame("Data"),
                        new JInternalFrame("Raw"),
                }, "Image");
        assertEquals("Image (2)", title);

        title = UIUtils.getUniqueFrameTitle(
                new JInternalFrame[]{
                        new JInternalFrame("Image (3)"),
                        new JInternalFrame("Data"),
                        new JInternalFrame("Raw"),
                }, "Image");
        assertEquals("Image", title);

        title = UIUtils.getUniqueFrameTitle(
                new JInternalFrame[]{
                        new JInternalFrame("Image"),
                        new JInternalFrame("Image (2)"),
                        new JInternalFrame("Image (3)"),
                }, "Image");
        assertEquals("Image (4)", title);

        title = UIUtils.getUniqueFrameTitle(
                new JInternalFrame[]{
                        new JInternalFrame("Image"),
                        new JInternalFrame("Image (2)"),
                        new JInternalFrame("Image (4)"),
                }, "Image");
        assertEquals("Image (3)", title);

        title = UIUtils.getUniqueFrameTitle(
                new JInternalFrame[]{
                        new JInternalFrame("Image"),
                        new JInternalFrame("Image (1)"),
                        new JInternalFrame("Image (2)"),
                        new JInternalFrame("Image (3)"),
                }, "Image");
        assertEquals("Image (4)", title);
    }

    public void testCreateSpinner_WithParameter() {
        final String labelname = "paramLabel";
        final ParamProperties properties = new ParamProperties(Integer.class, new Integer(3));
        properties.setLabel(labelname);
        final Parameter parameter = new Parameter("paramName", properties);

        final JSpinner spinner = UIUtils.createSpinner(parameter, new Integer(10), "#0");
        assertEquals(labelname, spinner.getName());
    }

    private void warnHeadless() {
        System.out.println("A " + UIUtilsTest.class + " test has not been performed: HeadlessException");
    }

}
