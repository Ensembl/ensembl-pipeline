package pathsteps.gui;

import javax.swing.*;
import pathsteps.model.*;

/**
 * This is the window which encloses the DataSourceConfigurationPanel
**/
public class LayoutConfigurationDialog extends JDialog {
  LayoutConfigurationPanel panel;
  
  public LayoutConfigurationDialog(PathStepsSwingView view){
    super(view.getRootFrame(), "Layout Parameters", false);
    panel = new LayoutConfigurationPanel(view);
    getContentPane().add(panel);
    view.connectToWindowCloseEventRouter(this,view.CANCEL_LAYOUT_PREFERENCES_KEY);
  }
  
  public void read(PathStepsModel model){
    getPanel().read(model);
  }
  
  public void update(PathStepsModel model){
    getPanel().update(model);
  }
  
  private LayoutConfigurationPanel getPanel(){
    return panel;
  }
}