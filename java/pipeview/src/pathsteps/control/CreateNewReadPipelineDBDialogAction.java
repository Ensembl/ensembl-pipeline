package pathsteps.control;

import java.util.*;
import java.sql.*;

import pathsteps.common.*;
import pathsteps.model.*;
import pathsteps.gui.*;

/**
 * Extend the model to include the new pipeline db info, filling in with
 * defaults/history as necessary.
 * Message the view to create a new pipeline db.
 * Exit - the view's information will be implicitly read from the model by the application.
**/
public class CreateNewReadPipelineDBDialogAction extends AAction{
  public CreateNewReadPipelineDBDialogAction(AEventHandler eventHandler){
    super(eventHandler);
  }
  
  protected void doAction(AView view, AModel aModel){
    Application application = view.getApplication();
    PathStepsModel model = (PathStepsModel)aModel;
    ModelElement dbDialogElement = null;
    ModelElement layoutDialogElement = null;
    String host;
    String port;
    String user;
    
    if(view.isReadPipelineDBDialogOpen()){
      if(isLoggingMedium()){
        logMedium("Dialog already open - not opening new dialog");
      }
      view.requestFocusOnPipelineDBDialog();
      return;
    }
    
    application.readHistory();

    if(model.getRootElement().getChildElement(PathStepsModel.READ_DB_DIALOG) == null){
      model.getRootElement().createChildElement(PathStepsModel.READ_DB_DIALOG);
      
      dbDialogElement = model.getRootElement().getChildElement(PathStepsModel.READ_DB_DIALOG);

      host = application.getHistory().getProperty(model.READ_DB_DIALOG_HOST);
      port = application.getHistory().getProperty(model.READ_DB_DIALOG_PORT);
      user = application.getHistory().getProperty(model.READ_DB_DIALOG_USER);

      if(host == null){
        host = application.getConfiguration().getProperty(model.READ_DB_DIALOG_HOST);
      }

      if(port == null){
        port = application.getConfiguration().getProperty(model.READ_DB_DIALOG_PORT);
      }

      if(user == null){
        user = application.getConfiguration().getProperty(model.READ_DB_DIALOG_USER);
      }

      dbDialogElement.addProperty(PathStepsModel.READ_DB_DIALOG_HOST, host);
      dbDialogElement.addProperty(PathStepsModel.READ_DB_DIALOG_PORT, port);
      dbDialogElement.addProperty(PathStepsModel.READ_DB_DIALOG_USER, user);
      dbDialogElement.addProperty(PathStepsModel.READ_DB_DIALOG_PASSWORD, "");
    }
    
    if(model.getRootElement().getChildElement(model.LAYOUT_DIALOG) == null){
      model.getRootElement().createChildElement(model.LAYOUT_DIALOG);
      layoutDialogElement = model.getRootElement().getChildElement(PathStepsModel.LAYOUT_DIALOG);
      layoutDialogElement.addProperty(model.LAYOUT_DIALOG_HORIZONTAL_SPACING, getStringFromHistoryOrConfig(model.LAYOUT_DIALOG_HORIZONTAL_SPACING));
      layoutDialogElement.addProperty(model.LAYOUT_DIALOG_VERTICAL_SPACING, getStringFromHistoryOrConfig(model.LAYOUT_DIALOG_VERTICAL_SPACING));
      layoutDialogElement.addProperty(model.LAYOUT_DIALOG_REPULSION_MULTIPLIER, getStringFromHistoryOrConfig(model.LAYOUT_DIALOG_REPULSION_MULTIPLIER));
      layoutDialogElement.addProperty(model.LAYOUT_DIALOG_ITERATES, getStringFromHistoryOrConfig(model.LAYOUT_DIALOG_ITERATES));
      layoutDialogElement.addProperty(model.LAYOUT_DIALOG_MOVEMENT_LIMIT, getStringFromHistoryOrConfig(model.LAYOUT_DIALOG_MOVEMENT_LIMIT));
      layoutDialogElement.addProperty(model.LAYOUT_DIALOG_GRAVITY, getStringFromHistoryOrConfig(model.LAYOUT_DIALOG_GRAVITY));
      layoutDialogElement.addProperty(model.LAYOUT_DIALOG_FIX_ROOTS, getStringFromHistoryOrConfig(model.LAYOUT_DIALOG_FIX_ROOTS));
      layoutDialogElement.addProperty(model.LAYOUT_DIALOG_SPRING_NATURAL_LENGTH, getStringFromHistoryOrConfig(model.LAYOUT_DIALOG_SPRING_NATURAL_LENGTH));
    }
    
    //
    //If we can find the db's, then populate them. Otherwise leave the list empty.
    createDatabaseList(model);

    if(isLoggingMedium()){
      logMedium("Opening pipeline db dialog");
    }
    
    view.openReadPipelineDBDialog();
    
    if(isLoggingMedium()){
      logMedium("Opened pipeline db dialog");
    }
    
    if(isLoggingLow()){
      logLow("Finished action");
    }
  }
  
  private void createDatabaseList(PathStepsModel model){
    
    if(isLoggingMedium()){
      logMedium("Fetching list of databases");
    }
    
    ModelElement element = model.getRootElement().getChildElement(PathStepsModel.READ_DB_DIALOG);
    ModelElement listElement = element.createChildElement(PathStepsModel.READ_DB_DIALOG_PIPELINE_DB_LIST);
    ArrayList orderedListOfNames = new ArrayList();
    listElement.addProperty(PathStepsModel.READ_DB_DIALOG_PIPELINE_DB_LIST, orderedListOfNames);
    if(isLoggingMedium()){
      logMedium("Adding "+orderedListOfNames+" to model");
    }

    String host = (String)element.getProperty(PathStepsModel.READ_DB_DIALOG_HOST);
    String port = (String)element.getProperty(PathStepsModel.READ_DB_DIALOG_PORT);
    String url = "jdbc:mysql://" + host + ":" + port + "/";
    String user = (String)element.getProperty(PathStepsModel.READ_DB_DIALOG_USER);
    String password = (String)element.getProperty(PathStepsModel.READ_DB_DIALOG_PASSWORD);
    
    String jdbcDriver = "org.gjt.mm.mysql.Driver";
    int hostIndex = -1;

    if(
      jdbcDriver == null ||
      jdbcDriver.trim().length() <= 0 ||
      host == null ||
      host.trim().length() <= 0 ||
      port == null ||
      port.trim().length() <= 0 ||
      user == null ||
      user.trim().length() <= 0
    ){
      
      return;
    }

    if(password == null){
      password = "";
    }
    
    try{
      Class.forName(jdbcDriver).newInstance();
    }catch(IllegalAccessException exception){
      
      if(isLoggingLow()){
        logLow("Could not access JDBC Driver constructor", exception);
      }
      
      model.getRootElement().addProperty(PathStepsModel.MESSAGE, "Could not access database with prior information");
      return;
    }catch(InstantiationException exception){

      if(isLoggingLow()){
        logLow("Could not create JDBC Driver", exception);
      }
      
      model.getRootElement().addProperty(PathStepsModel.MESSAGE, "Could not access database with prior information");
      return;
    }catch(ClassNotFoundException exception){

      if(isLoggingLow()){
        logLow("Could not access JDBC Driver class", exception);
      }
      
      model.getRootElement().addProperty(PathStepsModel.MESSAGE, "Could not access database with prior information");
      return;
    }//end try

    try{
      java.sql.Connection conn = DriverManager.getConnection(url,user,password);
      Statement statement = conn.createStatement();
      DatabaseMetaData metadata = conn.getMetaData();
      
      ResultSet databases = statement.executeQuery("show databases");
      String databaseName;
      
      while(databases.next()){
        databaseName = databases.getString(1);
        if(isLoggingHigh()){
          logHigh("fetched database name:"+databaseName);
        }
        orderedListOfNames.add(databaseName);
      }//end while
      Collections.sort(orderedListOfNames);
      
    }catch(SQLException exception){
      if(isLoggingLow()){
        logLow("Could not access JDBC Driver class", exception);
      }
      
      model.getRootElement().addProperty(PathStepsModel.MESSAGE, "SQL Problems reading pipline database");
      //throw new NonFatalAException("Problems reading pipeline database");
    }
    
    if(isLoggingMedium()){
      logMedium("Successfully finished fetching list of databases");
    }
  }
}
