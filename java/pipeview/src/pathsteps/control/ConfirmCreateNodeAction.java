package pathsteps.control;

import java.util.*;
import java.sql.*;

import pathsteps.common.*;
import pathsteps.model.*;
import pathsteps.gui.*;

/**
 * 1. Work out what the user actually wants to create </br>
 * 2. Do the db actions to insert a new analysis and rule_goal
 * 3. redo the graph-model so the new node will actually appear somewhere...
**/
public class ConfirmCreateNodeAction extends AAction{
  public ConfirmCreateNodeAction(AEventHandler eventHandler){
    super(eventHandler);
  }
  
  protected void doAction(AView view, AModel aModel){
    Application application = view.getApplication();
    PathStepsModel model = (PathStepsModel)aModel;
    ModelElement dialogElement = null;
    String name;
    String db;
    String dbVersion;
    String dbFile;
    String program;
    String programVersion;
    String programFile;
    String parameters;
    String module;
    String moduleVersion;
    String gffSource;
    String gffFeature;
    String type;
    
    dialogElement = model.getRootElement().getChildElement(PathStepsModel.CREATE_NODE_DIALOG);
    name = (String)dialogElement.getProperty(model.CREATE_NODE_DIALOG_NAME);
    db = (String)dialogElement.getProperty(model.CREATE_NODE_DIALOG_DB_FILE);
    dbFile = (String)dialogElement.getProperty(model.CREATE_NODE_DIALOG_DB_FILE);
    program = (String)dialogElement.getProperty(model.CREATE_NODE_DIALOG_PROGRAM);
    programVersion = (String)dialogElement.getProperty(model.CREATE_NODE_DIALOG_PROGRAM_VERSION);
    programFile = (String)dialogElement.getProperty(model.CREATE_NODE_DIALOG_PROGRAM_FILE);
    parameters = (String)dialogElement.getProperty(model.CREATE_NODE_DIALOG_PARAMETERS);
    module = (String)dialogElement.getProperty(model.CREATE_NODE_DIALOG_MODULE);
    moduleVersion = (String)dialogElement.getProperty(model.CREATE_NODE_DIALOG_MODULE_VERSION);
    gffSource = (String)dialogElement.getProperty(model.CREATE_NODE_DIALOG_GFF_SOURCE);
    gffFeature = (String)dialogElement.getProperty(model.CREATE_NODE_DIALOG_GFF_FEATURE);
    type = (String)dialogElement.getProperty(model.CREATE_NODE_DIALOG_INPUT_ID_TYPE);
    
    if(isLoggingMedium()){
      logMedium("creating a new node with name: "+name+" and type: "+type);
    }
    
    int analysisId = 
      createAnalysis(
        name,
        db,
        dbFile,
        program,
        programVersion,
        programFile,
        parameters,
        module,
        moduleVersion,
        gffSource,
        gffFeature
      );
    
    createInputIdTypeAnalysis(analysisId, type);
    
    createRuleGoalForAnalysis(name);
    
    //
    //Now forward to the read-db action to recreate the new model with the added element.
    if(isLoggingLow()){
      logLow("Forwarding - RereadSpecifiedDBAction");
    }
    
    getEventHandler().doActionForKey(AView.REFRESH_BUTTON_KEY);
    
    if(isLoggingLow()){
      logLow("Returning from forward - RereadSpecifiedDBAction");
    }
    
    if(isLoggingLow()){
      logLow("Finished ConfirmCreateNodeAction");
    }
  }
  
  private int createAnalysis(
    String name,
    String db,
    String dbFile,
    String program,
    String programVersion,
    String programFile,
    String parameters,
    String module,
    String moduleVersion,
    String gffSource,
    String gffFeature
  ){
    Connection connection = createNewConnectionFromSelectedDatabase((PathStepsModel)getModel());
    int createdId = 0;
    PreparedStatement statement;
    ResultSet insertSet;
    
    String sql = 
      "insert into analysis "+
      "set analysis_id = ?, "+
      "name = ?, "+
      "db = ?, "+
      "db_file = ?, "+
      "program = ?, "+
      "program_version = ?, "+
      "program_file = ?, "+
      "parameters = ?, "+
      "module = ?, "+
      "module_version = ?, "+
      "gff_source = ?, "+
      "gff_feature = ?";
    

    try{
      statement = connection.prepareStatement(sql);
      statement.setObject(1, null);
      statement.setString(2, db);
      statement.setString(3, dbFile);
      statement.setString(4, program);
      statement.setString(5, programVersion);
      statement.setString(6, parameters);
      statement.setString(7, module);
      statement.setString(8, moduleVersion);
      statement.setString(9, gffSource);
      statement.setString(10, gffFeature);

      int rowsInserted = statement.executeUpdate();
      rowsInserted = statement.executeUpdate();
      if(rowsInserted != 1){
        throw new NonFatalAException("Problem inserting analysis - statement executed but failed to insert");
      }
      
      insertSet = connection.createStatement().executeQuery("select last insert_id()");
      insertSet.next();
      createdId = insertSet.getInt(1);
      if(createdId <= 0){
        throw new NonFatalAException("Problem inserting analysis - insert created unnaceptable internal id");
      }

      connection.close();
    
    }catch(SQLException exception){
      throw new NonFatalAException("Problems inserting analysis: "+exception.getMessage(), exception);
    }
    
    return createdId;
  }
  
  private void createInputIdTypeAnalysis(int analysisId, String type){
    Connection connection = createNewConnectionFromSelectedDatabase((PathStepsModel)getModel());
    PreparedStatement statement;
    
    String sql = 
      "insert into input_id_type_analysis "+
      "set analysis_id = ?, "+
      "type = ?";

    try{
      statement = connection.prepareStatement(sql);
      statement.setInt(1, analysisId);
      statement.setString(2, type);

      int rowsInserted = statement.executeUpdate();
      rowsInserted = statement.executeUpdate();
      if(rowsInserted != 1){
        throw new NonFatalAException("Problem inserting input_id_type_analysis - statement executed but failed to insert");
      }
      
      connection.close();
      
    }catch(SQLException exception){
      throw new NonFatalAException("Problems inserting input_id_type_analysis: "+exception.getMessage(), exception);
    }
    
  }
  
  private void createRuleGoalForAnalysis(String analysisName){
    Connection connection = createNewConnectionFromSelectedDatabase((PathStepsModel)getModel());
    PreparedStatement statement;
    
    String sql = 
      "insert into rule_goal "+
      "set rule_id = ?, "+
      "goal = ?";

    try{
      statement = connection.prepareStatement(sql);
      statement.setObject(1, null);
      statement.setString(2, analysisName);

      int rowsInserted = statement.executeUpdate();
      rowsInserted = statement.executeUpdate();
      if(rowsInserted != 1){
        throw new NonFatalAException("Problem inserting rule_goal - statement executed but failed to insert");
      }

      connection.close();

    }catch(SQLException exception){
      throw new NonFatalAException("Problems inserting rule_goal: "+exception.getMessage(), exception);
    }
    
  }
}
