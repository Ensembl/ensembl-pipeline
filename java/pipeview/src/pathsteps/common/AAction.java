package pathsteps.common;
import pathsteps.model.*;

public abstract class AAction{
  private AView view;
  private AModel model;
  private AEventHandler eventHandler;
  
  public AAction(AEventHandler eventHandler){
    setView(eventHandler.getApplication().getView());
    setModel(eventHandler.getApplication().getModel());
    setEventHandler(eventHandler);
    
    if(getView()== null){
      throw new FatalAException("Attempt to create AAction with null AView");
    }
    
    if(getModel()== null){
      throw new FatalAException("Attempt to create AAction with null AModel");
    }
    
    if(getLogger()== null){
      throw new FatalAException("Attempt to create AAction with null ALogger");
    }
  }
    
  public void processAction(){
    if(getLogger().isLoggingLow()){
      getLogger().logLow("Processing Action "+getClass().getName());
    }
    
    doAction(getView(), getModel());
    
    if(getLogger().isLoggingLow()){
      getLogger().logLow("Finished Processing Action "+getClass().getName());
    }
  }
  
  protected AView getView(){
    return view;
  }
  
  protected void setView(AView view){
    this.view = view;
  }
  
  protected AModel getModel(){
    return model;
  }
  
  protected void setModel(AModel model){
    this.model = model;
  }
  
  public AEventHandler getEventHandler(){
    return eventHandler;
  }
  
  public void setEventHandler(AEventHandler newValue){
    eventHandler = newValue;
  }
  
  protected abstract void doAction(AView view, AModel model);
  
  protected String getStringFromHistoryOrConfig(String key){
    String aString = 
      getView()
        .getApplication()
        .getHistory()
        .getProperty(key);

    if(aString == null || aString.trim().length() <=0){
      aString = 
        getView()
          .getApplication()
          .getConfiguration()
          .getProperty(key);
    }
    
    if(aString == null || aString.trim().length() <=0){
      throw new FatalAException("You must provide a value in the application config for: "+key);
    }
    
    return aString;
  }
  
  protected int getNumberFromHistoryOrConfig(String key){    
    int number = 0;
    String numberString = getStringFromHistoryOrConfig(key);
    
    try{
      number = Integer.valueOf(numberString).intValue();
    }catch(NumberFormatException exception){
      throw new FatalAException(
        "The layout config item: "+key+" is "+numberString+
        " which is not a valid integer", 
        exception
      );
    }
    
    return number;
  }  
  
  protected java.sql.Connection createNewConnectionFromSelectedDatabase(pathsteps.model.PathStepsModel model){
    if(getLogger().isLoggingMedium()){
      getLogger().logMedium("Creating connection from database in READ_DB_DIALOG model");
    }
    
    java.sql.Connection connection;
    ModelElement element = model.getRootElement().getChildElement(PathStepsModel.READ_DB_DIALOG);

    if(element == null){
      throw new NonFatalAException("Attempt to access a database before initialising");
    }
    
    String host = (String)element.getProperty(PathStepsModel.READ_DB_DIALOG_HOST);
    String port = (String)element.getProperty(PathStepsModel.READ_DB_DIALOG_PORT);
    String database = (String)element.getProperty(PathStepsModel.READ_DB_DIALOG_PIPELINE_DB);
    String url = "jdbc:mysql://" + host + ":" + port + "/"+database;
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
      throw new NonFatalAException("DB Dialog model had insufficient information to create a db connection");
    }

    if(password == null){
      password = "";
    }
    
    try{
      
      Class.forName(jdbcDriver).newInstance();
      
    }catch(IllegalAccessException exception){
      
      if(getLogger().isLoggingLow()){
        getLogger().logLow("Could not access JDBC Driver constructor", exception);
      }
      throw new NonFatalAException("Could not create database connection: Could not access JDBC Driver constructor", exception);
      
    }catch(InstantiationException exception){
      if(getLogger().isLoggingLow()){
        getLogger().logLow("Could not create JDBC Driver", exception);
      }
      throw new NonFatalAException("Could not create database connection: Could not create JDBC Driver", exception);
      
    }catch(ClassNotFoundException exception){
      if(getLogger().isLoggingLow()){
        getLogger().logLow("Could not access JDBC Driver class", exception);
      }
      
      throw new NonFatalAException("Could not create database connection: Could not access JDBC Driver class", exception);
    }//end try

    try{
      
      connection = java.sql.DriverManager.getConnection(url,user,password);
      
    }catch(java.sql.SQLException exception){
      if(getLogger().isLoggingLow()){
        getLogger().logLow("Could not access JDBC Driver class", exception);
      }
      
      throw new NonFatalAException("Could not create database connection: "+exception.getMessage(),exception);
    }
    
    if(getLogger().isLoggingMedium()){
      getLogger().logMedium("Successfully created connection");
    }
    
    return connection;
  }

  public ALogger getLogger(){
    return getEventHandler().getLogger();
  }
  
  public boolean isLoggingLow(){
    return getLogger().isLoggingLow();
  }
  
  public void logLow(String message){
    getLogger().logLow(message);
  }
  
  public void logLow(String message, Throwable exception) {
    getLogger().logLow(message, exception);
  }
  
  public boolean isLoggingMedium() {
    return getLogger().isLoggingMedium();
  }
  
  public void logMedium(String message) {
    getLogger().logMedium(message);
  }
  
  public void logMedium(String message, Throwable exception) {
    getLogger().logMedium(message, exception);
  }
  
  public boolean isLoggingHigh(){
    return getLogger().isLoggingHigh();
  }
  
  public void logHigh(String message) {
    getLogger().logHigh(message);
  }
  
  public void logHigh(String message, Throwable exception) {
    getLogger().logHigh(message, exception);
  }  
}
