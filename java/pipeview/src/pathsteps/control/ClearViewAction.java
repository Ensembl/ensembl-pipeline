package pathsteps.control;
import pathsteps.common.*;
import pathsteps.model.*;
import pathsteps.gui.*;

/**
 * When the user types into the db-change field, the databases in the dropdown
 * should be blanked out: this is what this action does.
**/
public class ClearViewAction extends AAction{
  public ClearViewAction(AEventHandler eventHandler){
    super(eventHandler);
  }
  
  protected void doAction(AView view, AModel aModel){
    view.reshow();
  }
}
