package pathsteps.common;
import java.awt.event.*;
import java.util.*;
import java.util.logging.*;
import java.io.*;

/**
 * The path steps application: has the main method, and holds the communication between
 * all components
**/
public class Application{
    
  public static String PIPELINE_TOOL_ROOT = "PIPELINE_TOOL_ROOT";
  public static String DEFAULT_CONFIGURATION_FILE_NAME = "pathsteps.conf";
  public static String DEFAULT_GRAPH_LAYOUT_CONFIGURATION_FILE_NAME = "graphlayout.conf";
  public static String DEFAULT_HISTORY_FILE_NAME = "pathsteps.history";
  public static String LOGGING_CONFIG_FILE_NAME = "logging.properties";
  public static String JAVA_LOGGING_PROPERTIES_FILE_NAME = "java.util.logging.config.file";
  public static String DEFAULT_TEST_CONFIGURATION_FILE_NAME  = "test.conf";
  
  private AView view;
  private AModel model;
  
  private AEventHandler eventHandler;
  
  private ALogger eventLogger;
  private ALogger applicationLogger;
  private ALogger viewLogger;
  private ALogger modelLogger;
  
  private String applicationRoot;
  
  private Properties graphLayoutConfiguration = new Properties();
  private Properties configuration = new Properties();
  private Properties history = new Properties();
  
  private boolean testMode = false;
  
  public static void main(String[] args){
    Application application = null;
    String root = System.getProperty("PIPELINE_TOOL_ROOT");

    try{
      if(root == null){
        throw new FatalAException("You must set a system variable PIPELINE_TOOL_ROOT");
      }
      
      //The app is just an empty shell - a backbone hanging onto everything else.
      application = new Application();
      
      //make sure that the system property is set.
      application.setApplicationRoot(root);
      
      //if we have a logging.properties file in the system root, then we can configure
      //a java logger from it. Otherwise use a simple console logger.
      application.createLoggers();
      
      //get the config file
      application.readConfiguration();
      
      //See if we've embedded a test configuration file.
      application.readTestConfiguration();
      
      
      application.readHistory();
      
      if(!application.isTestMode()){
        application.setView(new pathsteps.gui.PathStepsSwingView());
      }else{
        application.setView(new pathsteps.test.PathStepsTestView());
      }
      
      application.setModel(new pathsteps.model.PathStepsModel());
      
      application.setEventHandler(new pathsteps.control.PathStepsEventHandler());
      application.getEventHandler().setApplication(application);
      application.getEventHandler().setLogger(application.getEventLogger());
      application.getEventHandler().initialise();
      
      application.getView().setApplication(application);
      application.getView().setLogger(application.getViewLogger());
      application.getView().initialise();
      application.getView().show();
      
      if(application.isTestMode()){
        application.runTests();
      }
      
    }catch(Throwable exception){
      //
      //Have to make a raw reference to swing code here: can't use the view if the view is
      //refusing to be created. This should be the ONLY place in the app that this happens.
      exception.printStackTrace();
      javax.swing.JOptionPane.showMessageDialog(null, "Fatal Error During Startup: \n"+exception.getMessage());
    }//end try
  }

  /**
   * If there is a logging.properties file in the application root, then use it to configure
   * a set of java-standard loggers. Otherwise, use 'SimpleLogger's which just write prints
   * to the console.
  **/
  private void createLoggers(){
    //
    //This weird little line has the effect of LOADING the ALogLevel class, 
    //which instantiates the array of 'allowed' LogLevels inside LogLevel, which
    //in turn means that the LogManager can happily set itself to a log level
    //of LOW/MEDIUM/HIGH, instead of one of the standard log levels...
    ALogLevel level = new ALogLevel("Hello", 1);
    
    String fileName = getLoggerConfigurationFileName();
    File configurationFile;
    LogManager manager;
    
    configurationFile = new File(fileName);
    if(!configurationFile.exists()){
      setEventLogger(new SimpleLogger(ALogRecord.CONTROL));
      setApplicationLogger(new SimpleLogger(ALogRecord.APP));
      setViewLogger(new SimpleLogger(ALogRecord.VIEW));
    }else{
      
      System.setProperty(JAVA_LOGGING_PROPERTIES_FILE_NAME, fileName);
      manager = LogManager.getLogManager();
      try{
        manager
          .readConfiguration(
            new FileInputStream(fileName)
          );
      }catch(IOException exception){
        throw new FatalAException("Could not successfully open the logging properties file: "+fileName);
      }
      
      setEventLogger(new JavaLogger(ALogRecord.CONTROL));
      setApplicationLogger(new JavaLogger(ALogRecord.APP));
      setViewLogger(new JavaLogger(ALogRecord.VIEW));
      setModelLogger(new JavaLogger(ALogRecord.MODEL));
    }
  }

  private void setTestMode(boolean newValue){
    testMode = newValue;
  }

  private boolean isTestMode(){
    return testMode;
  }
  
  public AView getView(){
    return view;
  }
  
  private void setView(AView view){
    this.view = view;
  }
  
  public AModel getModel(){
    return this.model;
  }
  
  private void setModel(AModel model){
    this.model = model;
  }
  
  private AEventHandler getEventHandler(){
    return eventHandler;
  }
  
  private void setEventHandler(AEventHandler eventHandler){
    this.eventHandler = eventHandler;
  }
  
  private ALogger getEventLogger(){
    return eventLogger;
  }
  
  private void setEventLogger(ALogger logger){
    eventLogger = logger;
  }
  
  private ALogger getApplicationLogger(){
    return applicationLogger;
  }
  
  private void setApplicationLogger(ALogger logger){
    applicationLogger = logger;
  }
  
  private ALogger getViewLogger(){
    return viewLogger;
  }

  private void setViewLogger(ALogger logger){
    viewLogger = logger;
  }
  
  private ALogger getModelLogger(){
    return modelLogger;
  }

  private void setModelLogger(ALogger logger){
    modelLogger = logger;
  }
  
  /**
   * This is the signal that passes between the EventRouters in the gui-view
   * and the EventHandler's actions
  **/
  public void notifyEventForKey(String key){
    getEventHandler().notifyEventForKey(key);
  }

  private void readTestConfiguration(){
    String fileName = getTestConfigurationFileName();
    File configurationFile = new File(fileName);
    if(!configurationFile.exists()){
      setTestMode(false);
    }else{
      setTestMode(true);
    }
  }
  
  private void readConfiguration(){
    String fileName = getConfigurationFileName();
    File configurationFile;
    
    try{
      
      configurationFile = new File(fileName);
      if(!configurationFile.exists()){
        throw new FatalAException("You must have an application configuration file: "+fileName);
      }
      getConfiguration().load(new FileInputStream(fileName));
      
    }catch(IOException exception){
      throw new FatalAException("Could not read application properties file: "+fileName, exception);
    }
  }
  
  public void readHistory(){
    String fileName = getHistoryFileName();
    try{
      getHistory().load(new FileInputStream(fileName));
    }catch(IOException exception){
      throw new FatalAException("Could not read application properties file: "+fileName, exception);
    }
  }
  
  public void writeHistory(Properties history){
    String fileName = getHistoryFileName();
    if(getApplicationLogger().isLoggingLow()){
      getApplicationLogger().logLow("writing history to file :"+fileName+" contents: "+history);
    }
    try{
      history.store(new FileOutputStream(fileName), "Application History");
    }catch(IOException exception){
      throw new FatalAException("Could not write application history file: "+fileName, exception);
    }
  }
  
  public void writeGraphLayoutConfiguration(Properties history){
    String fileName = getGraphLayoutConfigurationFileName();
    if(getApplicationLogger().isLoggingLow()){
      getApplicationLogger().logLow("writing graph layout config to file :"+fileName+" contents: "+history);
    }
    try{
      history.store(new FileOutputStream(fileName), "Graph Layout Configuration");
    }catch(IOException exception){
      throw new FatalAException("Could not graph layout configuration for file: "+fileName, exception);
    }
  }
  
  public Properties readGraphLayoutConfiguration(){
    String fileName = getGraphLayoutConfigurationFileName();
    File configFile = new File(fileName);
    
    if(configFile.exists()){
      if(getApplicationLogger().isLoggingLow()){
        getApplicationLogger().logLow("Reading graph layout config from file :"+fileName);
      }

      try{
        getGraphLayoutConfiguration().load(new FileInputStream(fileName));
      }catch(IOException exception){
        throw new FatalAException("Could not graph layout configuration for file: "+fileName, exception);
      }
    }
    return getGraphLayoutConfiguration();
  }
  
  private void setGraphLayoutConfiguration(Properties properties){
    graphLayoutConfiguration = properties;
  }
  
  public Properties getGraphLayoutConfiguration(){
    return graphLayoutConfiguration;
  }
  
  private void setConfiguration(Properties properties){
    configuration = properties;
  }
  
  private void setHistory(Properties properties){
    history = properties;
  }
  
  public Properties getConfiguration(){
    return configuration;
  }
  
  public Properties getHistory(){
    return history;
  }
  
  private String getConfigurationFileName(){
    return getApplicationRoot()+"/"+DEFAULT_CONFIGURATION_FILE_NAME;
  }
  
  private String getGraphLayoutConfigurationFileName(){
    return getApplicationRoot()+"/"+DEFAULT_GRAPH_LAYOUT_CONFIGURATION_FILE_NAME;
  }
  
  private String getTestConfigurationFileName(){
    return getApplicationRoot()+"/"+DEFAULT_TEST_CONFIGURATION_FILE_NAME;
  }
  
  private String getHistoryFileName(){
    return getApplicationRoot()+"/"+DEFAULT_HISTORY_FILE_NAME;
  }

  private String getLoggerConfigurationFileName(){
    return getApplicationRoot()+"/"+LOGGING_CONFIG_FILE_NAME;
  }
  
  private void setApplicationRoot(String root){
    applicationRoot = root;
  }
  
  public String getApplicationRoot(){
    return applicationRoot;
  }
  
  private void runTests(){
    pathsteps.test.TestRunner runner = 
      new pathsteps.test.TestRunner(
        this, 
        (pathsteps.test.PathStepsTestView)getView(), 
        (pathsteps.model.PathStepsModel)getModel()
      );
    
    runner.run();
  }

  public boolean isLoggingLow(){
    return getApplicationLogger().isLoggingLow();
  }
  
  public void logLow(String message){
    getApplicationLogger().logLow(message);
  }
  
  public void logLow(String message, Throwable exception) {
    getApplicationLogger().logLow(message, exception);
  }
  
  public boolean isLoggingMedium() {
    return getApplicationLogger().isLoggingMedium();
  }
  
  public void logMedium(String message) {
    getApplicationLogger().logMedium(message);
  }
  
  public void logMedium(String message, Throwable exception) {
    getApplicationLogger().logMedium(message, exception);
  }
  
  public boolean isLoggingHigh(){
    return getApplicationLogger().isLoggingHigh();
  }
  
  public void logHigh(String message) {
    getApplicationLogger().logHigh(message);
  }
  
  public void logHigh(String message, Throwable exception) {
    getApplicationLogger().logHigh(message, exception);
  }  
}