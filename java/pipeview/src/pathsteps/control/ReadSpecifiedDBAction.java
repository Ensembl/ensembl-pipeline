package pathsteps.control;

import pathsteps.common.*;
import pathsteps.model.*;

import java.util.*;
import java.sql.*;

/**
 * Read information from chosen database (communicated via the model)
 * and close the pipeline db dialog after we're done.
**/
public class ReadSpecifiedDBAction extends AbstractReadDBAction{

  public ReadSpecifiedDBAction(AEventHandler eventHandler){
    super(eventHandler);
  }
  
  protected void doAction(AView view, AModel aModel){
    super.doAction(view, aModel);
    
    view.closeReadPipelineDBDialog();
    
    if(isLoggingMedium()){
      logLow("Closed Read Pipeline DB Dialog");
    }
  }

}
