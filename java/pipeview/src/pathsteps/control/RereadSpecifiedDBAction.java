package pathsteps.control;

import pathsteps.common.*;
import pathsteps.model.*;

import java.util.*;
import java.sql.*;

/**
 * Reads the chosen db parameters, queries the db and populates
 * the results into the application model, but does NOT close the
 * pipeline db dialog - this is not open!
**/
public class RereadSpecifiedDBAction extends AbstractReadDBAction{

  public RereadSpecifiedDBAction(AEventHandler eventHandler){
    super(eventHandler);
  }
  
  protected void doAction(AView view, AModel aModel){
    super.doAction(view, aModel);
  }
  
}
