package pathsteps.control;
import pathsteps.common.*;
import pathsteps.model.*;
import java.sql.*;
import java.util.*;

public class FindAllDatabasesAction extends AAction{
  public FindAllDatabasesAction(AEventHandler eventHandler){
    super(eventHandler);
  }
  
  protected void doAction(AView view, AModel model){
    if(isLoggingMedium()){
      logLow("Started Find All DBases action");
    }
    
    createDatabaseList((PathStepsModel)model);
    
    if(isLoggingMedium()){
      logLow("Finished Find All DBases action");
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

      if(isLoggingMedium()){
        logMedium("Adding "+orderedListOfNames+" to model");
      }
    }catch(SQLException exception){
      if(isLoggingMedium()){
        logMedium("Could not access JDBC Driver class", exception);
      }
      
      model.getRootElement().addProperty(PathStepsModel.MESSAGE, "SQL Problems reading pipline database");
      return;
    }catch(Throwable exception){
      if(isLoggingMedium()){
        logMedium("Caught an unexpected exception when trying to connect to db...");
      }
      throw new NonFatalAException("Caught an unexpected exception when trying to connect to db:"+exception.getMessage(), exception);
    }
    
    if(isLoggingMedium()){
      logMedium("Successfully finished fetching list of databases");
    }
  }
  
}
