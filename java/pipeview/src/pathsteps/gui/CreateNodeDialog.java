package pathsteps.gui;

import javax.swing.*;
import pathsteps.model.*;

/**
 * This is the window which encloses the CreateNodePanel
**/
public class CreateNodeDialog extends JDialog {
  CreateNodePanel panel;
  
  public CreateNodeDialog(PathStepsSwingView view){
    super(view.getRootFrame(), "Create New Node", false);
    panel = new CreateNodePanel(view);
    getContentPane().add(panel);
    view.connectToWindowCloseEventRouter(this,view.CANCEL_NODE_DIALOG_KEY);
  }
  
  public void read(PathStepsModel model){
    getPanel().read(model);
  }
  
  public void update(PathStepsModel model){
    getPanel().update(model);
  }
  
  private CreateNodePanel getPanel(){
    return panel;
  }
}