package pathsteps.control;
import pathsteps.common.*;
import pathsteps.model.*;

/**
 * When the user types into the db-change field, the databases in the dropdown
 * should be blanked out: this is what this action does.
**/
public class CancelLayoutConfigurationAction extends AAction{
  public CancelLayoutConfigurationAction(AEventHandler eventHandler){
    super(eventHandler);
  }
  
  protected void doAction(AView view, AModel aModel){
    view.closeLayoutConfigurationDialog();
  }
}
