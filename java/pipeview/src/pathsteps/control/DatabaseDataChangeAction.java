package pathsteps.control;
import pathsteps.common.*;
import pathsteps.model.*;
import java.util.*;

/**
 * When the user types into the db-change field, the databases in the dropdown
 * should be blanked out: this is what this action does.
**/
public class DatabaseDataChangeAction extends AAction{
  public DatabaseDataChangeAction(AEventHandler eventHandler){
    super(eventHandler);
  }
  
  protected void doAction(AView view, AModel model){
    logLow("Starting database data changed action");
    PathStepsModel theModel = (PathStepsModel)model;
    ModelElement dialogModel = theModel.getRootElement().getChildElement(theModel.READ_DB_DIALOG);
    ModelElement child = dialogModel.getChildElement(theModel.READ_DB_DIALOG_PIPELINE_DB_LIST);
    child.addProperty(theModel.READ_DB_DIALOG_PIPELINE_DB_LIST, new ArrayList());
    logLow("Finished database data changed action");
  }
}
