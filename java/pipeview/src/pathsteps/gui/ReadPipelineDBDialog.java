package pathsteps.gui;

import javax.swing.*;
import pathsteps.model.*;

/**
 * This is the window which encloses the DataSourceConfigurationPanel
**/
public class 
  ReadPipelineDBDialog
extends 
  JDialog
{
  DataSourceConfigurationPanel panel;
  
  public ReadPipelineDBDialog(PathStepsSwingView view){
    super(view.getRootFrame(), "Find databases", false);
    panel = new DataSourceConfigurationPanel(view);
    getContentPane().add(panel);
    view.connectToWindowCloseEventRouter(this,view.CANCEL_READ_KEY);
  }
  
  public void read(PathStepsModel model){
    getPanel().read(model);
  }
  
  public void update(PathStepsModel model){
    getPanel().update(model);
  }
  
  private DataSourceConfigurationPanel getPanel(){
    return panel;
  }
}