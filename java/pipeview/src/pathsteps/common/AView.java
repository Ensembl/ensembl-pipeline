package pathsteps.common;

/**
 * The interface from the application to its view. Implementations can be
 * graphical or test-based. The separation via the interface suggests a useful
 * test method for the gui application as a whole.
**/
public interface AView{
  public static String SAVE_BUTTON_KEY = "SAVE_BUTTON_KEY";
  public static String RESET_BUTTON_KEY = "RESET_BUTTON_KEY";
  public static String REFRESH_BUTTON_KEY = "CLEAR_BUTTON_KEY";
  public static String LAYOUT_BUTTON_KEY = "LAYOUT_BUTTON_KEY";
  public static String SET_NEW_PIPELINE_DB_KEY = "SET_NEW_PIPELINE_DB_KEY";
  public static String SET_NEW_LAYOUT_PREFERENCES_KEY = "SET_NEW_LAYOUT_PREFERENCES_KEY";
  public static String READ_NEW_PIPELINE_DB_KEY = "READ_NEW_PIPELINE_DB_KEY";
  public static String CLOSE_ROOT_FRAME_KEY = "CLOSE_ROOT_FRAME_KEY";
  public static String DATABASE_DATA_CHANGE_KEY = "DATABASE_DATA_CHANGE_KEY";
  public static String FIND_ALL_DATABASES_KEY = "FIND_ALL_DATABASES_KEY";
  public static String READ_SPECIFIED_DB_KEY = "READ_SPECIFIED_DB_KEY";
  public static String CANCEL_READ_KEY = "CANCEL_READ_KEY";
  public static String SET_LAYOUT_PREFERENCES_KEY = "SET_LAYOUT_PREFERENCES_KEY";
  public static String CANCEL_LAYOUT_PREFERENCES_KEY = "CANCEL_LAYOUT_PREFERENCES_KEY";
  public static String SELECT_GRAPH_NODE_PARENTS_KEY  = "SELECT_GRAPH_NODE_KEY";
  
  public static String CREATE_NODE_DIALOG_KEY  = "CREATE_NODE_DIALOG_KEY";
  public static String UPDATE_NODE_DIALOG_KEY  = "UPDATE_NODE_DIALOG_KEY";
  public static String DELETE_NODE_DIALOG_KEY  = "DELETE_NODE_DIALOG_KEY";
  public static String CONNECT_NODE_DIALOG_KEY  = "CONNECT_NODE_DIALOG_KEY";
  
  public static String CONFIRM_CREATE_NODE_KEY  = "CONFIRM_CREATE_NODE_KEY";
  public static String CONFIRM_UPDATE_NODE_KEY  = "CONFIRM_UPDATE_NODE_KEY";
  public static String CONFIRM_DELETE_NODE_KEY  = "CONFIRM_DELETE_NODE_KEY";
  public static String CONFIRM_CONNECT_NODE_KEY  = "CONFIRM_CONNECT_NODE_KEY";
  
  public static String SAVE_GRAPH_LAYOUT_CONFIGURATION_KEY  = "SAVE_GRAPH_LAYOUT_CONFIGURATION_KEY";
  
  public static String CANCEL_NODE_DIALOG_KEY  = "CANCEL_NODE_DIALOG_KEY";

  public void read(AModel model);
  public void update(AModel model);
  
  public Application getApplication();
  
  public void setApplication(Application application);
  
  public void initialise();
  
  public void show();

  public void showMessage(String message);
  
  public void showFatalDialog(String message);
  
  public void shutDown();
  
  public ALogger getLogger();
  public void setLogger(ALogger logger);
  
  public boolean isLoggingLow();
  public void logLow(String message);
  public void logLow(String message, Throwable exception);
  
  public boolean isLoggingMedium();
  public void logMedium(String message);
  public void logMedium(String message, Throwable exception);
  
  public boolean isLoggingHigh();
  public void logHigh(String message);
  public void logHigh(String message, Throwable exception);
  
  //app specific things
  public boolean isReadPipelineDBDialogOpen();
  public void openReadPipelineDBDialog();
  public void closeReadPipelineDBDialog();
  public void requestFocusOnPipelineDBDialog();
  
  public boolean isLayoutConfigurationDialogOpen();
  public void openLayoutConfigurationDialog();
  public void closeLayoutConfigurationDialog();
  public void requestFocusOnLayoutConfigurationDialog();
  
  public boolean isCreateNodeDialogOpen();
  public void openCreateNodeDialog();
  public void closeCreateNodeDialog();
  public void requestFocusOnCreateNodeDialog();
  
  public String getNameOfSelectedNode();
  
  public void reshow();
  
  public void applyGraphLayout();
}//end AView
