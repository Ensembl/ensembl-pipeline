package pathsteps.control;

import pathsteps.common.*;
import pathsteps.model.*;

import java.util.*;
import java.sql.*;

/**
 * When the user types into the db-change field, the databases in the dropdown
 * should be blanked out: this is what this action does.
**/
public class ApplyGraphLayoutAction extends AAction{

  public ApplyGraphLayoutAction(AEventHandler eventHandler){
    super(eventHandler);
  }
  
  protected void doAction(AView view, AModel aModel){
    //
    //First write out the old graph layout to disk, so that the user can come back to it
    //if necessary.
    //getEventHandler().doActionForKey(view.SAVE_GRAPH_LAYOUT_CONFIGURATION_KEY);
    
    //
    //Now just prompt the view to do a relayout. 
    view.applyGraphLayout();
  }
}
