package com.bc.swing.dock;

import com.jidesoft.docking.DockContext;
import com.jidesoft.docking.DockableFrame;
import com.jidesoft.docking.DockingManager;

import javax.swing.AbstractAction;
import javax.swing.Icon;
import java.awt.Component;
import java.awt.Rectangle;
import java.awt.Window;
import java.awt.event.ActionEvent;
import java.awt.event.WindowListener;

/**
 *
 * @author Marco Peters
 * @version $Revision: 1.1 $ $Date: 2009-04-28 14:17:17 $
 * @since BEAM 4.5
 */
public class FloatingDockableFrame implements FloatingComponent {

    private static final int CONTENT_INDEX = 0;

    private DockableFrame dockableFrame;
    private DockableComponent originator;
    private final Window owner;
    private final DockingManager dockManager;

    public FloatingDockableFrame(Window owner, final DockingManager dockManager) {
        this.owner = owner;
        this.dockManager = dockManager;
        dockableFrame = new DockableFrame();
        dockableFrame.setFloatable(true);
        dockableFrame.setDockable(true);
        dockableFrame.getContext().setInitMode(DockContext.STATE_FLOATING);

        dockableFrame.setCloseAction(new AbstractAction() {

            public void actionPerformed(ActionEvent e) {
                if (getOriginator() != null) {
                    getOriginator().setDocked(true);
                }
            }
        });
        dockableFrame.setVisible(true);
    }

    public static FloatingComponentFactory getFactory(DockingManager dockManager) {
        return new Factory(dockManager);
    }

    public String getTitle() {
        return dockableFrame.getTitle();
    }

    public void setTitle(String title) {
        dockableFrame.setTitle(title);
    }

    public Icon getIcon() {
        return dockableFrame.getFrameIcon();
    }

    public void setIcon(Icon icon) {
        dockableFrame.setFrameIcon(icon);
    }

    public Component getContent() {
        return dockableFrame.getContentPane().getComponent(CONTENT_INDEX);
    }

    public void setContent(Component component) {
        dockableFrame.getContentPane().add(component, CONTENT_INDEX);
    }

    public DockableComponent getOriginator() {
        return originator;
    }

    public void setOriginator(DockableComponent originator) {
        this.originator = originator;
    }

    public Rectangle getBounds() {
        return dockableFrame.getBounds();
    }

    public void setBounds(Rectangle bounds) {
        dockableFrame.setBounds(bounds);
    }

    public void show() {
        dockableFrame.setKey(owner.getName() + "." + getTitle());
        dockManager.addFrame(dockableFrame);
        dockableFrame.setVisible(true);
    }

    public void close() {
        dockableFrame.setVisible(false);
        dockManager.removeFrame(dockableFrame.getKey());
    }

    public void addWindowListener(final WindowListener l) {

    }

    public void removeWindowListener(WindowListener l) {
    }

    private static class Factory implements FloatingComponentFactory {

        private DockingManager dockManager;

        private Factory(DockingManager dockManager) {
            this.dockManager = dockManager;
        }

        public FloatingComponent createFloatingComponent(Window owner) {

            return new FloatingDockableFrame(owner, dockManager);
        }
    }

}
