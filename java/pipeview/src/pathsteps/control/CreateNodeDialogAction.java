package pathsteps.control;

import java.util.*;
import java.sql.*;

import pathsteps.common.*;
import pathsteps.model.*;
import pathsteps.gui.*;

/**
 * 1. Extend the model to include the dialog model behind the gui.</br>
 * 2. Message the view to create a new "create node" dialog <br>
 * 3. Exit - the view's information will be implicitly read from the model by the application.
**/
public class CreateNodeDialogAction extends AAction{
  public CreateNodeDialogAction(AEventHandler eventHandler){
    super(eventHandler);
  }
  
  protected void doAction(AView view, AModel aModel){
    Application application = view.getApplication();
    PathStepsModel model = (PathStepsModel)aModel;
    ModelElement dialogElement = null;
    
    if(view.isCreateNodeDialogOpen()){
      if(isLoggingMedium()){
        logMedium("Create node dialog already open - not opening new dialog");
      }
      view.requestFocusOnCreateNodeDialog();
      return;
    }
          
    dialogElement = model.getRootElement().getChildElement(PathStepsModel.CREATE_NODE_DIALOG);
    
    if(dialogElement == null){
      if(isLoggingMedium()){
        logMedium("NO create node dialog element found - creating one");
      }
      
      model.getRootElement().createChildElement(PathStepsModel.CREATE_NODE_DIALOG);
      dialogElement = model.getRootElement().getChildElement(PathStepsModel.CREATE_NODE_DIALOG);
      
      application.readHistory(); //picks up the user's last choices from history instead of internal state.
    }else{
      if(isLoggingMedium()){
        logMedium("Create node dialog element found!");
      }
    }

    createInputIdTypeList(model);

    if(isLoggingMedium()){
      logMedium("populated input id list");
    }
    
    view.openCreateNodeDialog();
    
    if(isLoggingLow()){
      logLow("Finished action");
    }
  }
  
  private void createInputIdTypeList(PathStepsModel model){
    java.sql.Connection connection = createNewConnectionFromSelectedDatabase(model);
    ResultSet set;
    List inputIdTypeList = new ArrayList();
    ModelElement element = model.getRootElement().getChildElement(PathStepsModel.CREATE_NODE_DIALOG);
    String sql = "select distinct(input_id_type) from input_id_type_analysis";
    element.addProperty(PathStepsModel.CREATE_NODE_DIALOG_INPUT_ID_TYPE_LIST, inputIdTypeList);
    String name;
    
    if(isLoggingHigh()){
      logHigh("Running query to populate input id types list ");
    }
    
    try{
      
      set = 
        connection.createStatement().executeQuery(sql);
   
      if(isLoggingHigh()){
        logHigh("SQL: "+sql);
      }
      
      while(set.next()){
        name = set.getString(1);
        inputIdTypeList.add(name);
        if(isLoggingHigh()){
          logHigh("adding : "+name+" to allowed input id type list");
        }
      }
      
      Collections.sort(inputIdTypeList);
      if(isLoggingHigh()){
        logHigh("Size of input_id_type_list: "+inputIdTypeList.size());
      }
      
    }catch(SQLException exception){
      throw new NonFatalAException("Problems issueing SQL "+exception.getMessage(),exception);
    }
  }
}
