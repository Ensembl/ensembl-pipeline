package pathsteps.control;

import java.util.*;
import java.sql.*;

import pathsteps.common.*;
import pathsteps.model.*;
import pathsteps.gui.*;

/**
 * 1. Extend the model to include the layout configuration info, filling in with
 * defaults/history as necessary.</br>
 * 2. Message the view to create a new layout configuration dialog <br>
 * 3. Exit - the view's information will be implicitly read from the model by the application.
**/
public class CreateNewLayoutConfigurationDialogAction extends AAction{
  public CreateNewLayoutConfigurationDialogAction(AEventHandler eventHandler){
    super(eventHandler);
  }
  
  protected void doAction(AView view, AModel aModel){
    Application application = view.getApplication();
    PathStepsModel model = (PathStepsModel)aModel;
    ModelElement layoutDialogElement = null;

    if(view.isLayoutConfigurationDialogOpen()){
      if(isLoggingMedium()){
        logMedium("Layout configuration dialog already open - not opening new dialog");
      }
      view.requestFocusOnLayoutConfigurationDialog();
      return;
    }
          
    layoutDialogElement = model.getRootElement().getChildElement(PathStepsModel.LAYOUT_DIALOG);
    if(layoutDialogElement == null){
      if(isLoggingMedium()){
        logMedium("NO layout dialog element found - creating one");
      }
      
      model.getRootElement().createChildElement(PathStepsModel.LAYOUT_DIALOG);
      layoutDialogElement = model.getRootElement().getChildElement(PathStepsModel.LAYOUT_DIALOG);
      application.readHistory(); //picks up the user's last choices from history instead of internal state.
      
      layoutDialogElement.addProperty(model.LAYOUT_DIALOG_HORIZONTAL_SPACING, getStringFromHistoryOrConfig(model.LAYOUT_DIALOG_HORIZONTAL_SPACING));
      layoutDialogElement.addProperty(model.LAYOUT_DIALOG_VERTICAL_SPACING, getStringFromHistoryOrConfig(model.LAYOUT_DIALOG_VERTICAL_SPACING));
      layoutDialogElement.addProperty(model.LAYOUT_DIALOG_REPULSION_MULTIPLIER, getStringFromHistoryOrConfig(model.LAYOUT_DIALOG_REPULSION_MULTIPLIER));
      layoutDialogElement.addProperty(model.LAYOUT_DIALOG_ITERATES, getStringFromHistoryOrConfig(model.LAYOUT_DIALOG_ITERATES));
      layoutDialogElement.addProperty(model.LAYOUT_DIALOG_MOVEMENT_LIMIT, getStringFromHistoryOrConfig(model.LAYOUT_DIALOG_MOVEMENT_LIMIT));
      layoutDialogElement.addProperty(model.LAYOUT_DIALOG_GRAVITY, getStringFromHistoryOrConfig(model.LAYOUT_DIALOG_GRAVITY));
      layoutDialogElement.addProperty(model.LAYOUT_DIALOG_FIX_ROOTS, getStringFromHistoryOrConfig(model.LAYOUT_DIALOG_FIX_ROOTS));
      layoutDialogElement.addProperty(model.LAYOUT_DIALOG_SHOW_JOB_DETAIL, getStringFromHistoryOrConfig(model.LAYOUT_DIALOG_SHOW_JOB_DETAIL));
      layoutDialogElement.addProperty(model.LAYOUT_DIALOG_SPRING_NATURAL_LENGTH, getStringFromHistoryOrConfig(model.LAYOUT_DIALOG_SPRING_NATURAL_LENGTH));
    }else{
      if(isLoggingMedium()){
        logMedium("Layout dialog element found!");
      }
    }

    if(isLoggingMedium()){
      logMedium("Opening pipeline db dialog");
    }
    
    view.openLayoutConfigurationDialog();
    
    if(isLoggingMedium()){
      logMedium("Opened pipeline db dialog");
    }
  }
}
