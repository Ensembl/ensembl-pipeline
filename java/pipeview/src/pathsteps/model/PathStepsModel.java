package pathsteps.model;
import pathsteps.common.*;

public class PathStepsModel implements AModel{
  public static String ROOT = "ROOT";
  public static String PATH_STEPS_PANEL = "PATH_STEPS_PANEL";
  public static String PATH_STEPS_PANEL_SHOW_JOB_DETAIL = "PATH_STEPS_PANEL_SHOW_JOB_DETAIL";
  public static String PATH_STEPS_PANEL_GRAPH_LAYOUT_CONFIGURATION = "PATH_STEPS_PANEL_GRAPH_LAYOUT_CONFIGURATION";
  public static String PATH_STEPS_PANEL_ALL_NODES = "PATH_STEPS_PANEL_ALL_NODES";
  public static String PATH_STEPS_PANEL_ALL_NODES_MAP = "PATH_STEPS_PANEL_ALL_NODES_MAP";
  
  public static String MESSAGE = "MESSAGE";
  
  public static String READ_DB_DIALOG  = "READ_DB_DIALOG";
  public static String READ_DB_DIALOG_HOST = "READ_DB_DIALOG_HOST";
  public static String READ_DB_DIALOG_PORT = "READ_DB_DIALOG_PORT";
  public static String READ_DB_DIALOG_USER = "READ_DB_DIALOG_USER";
  public static String READ_DB_DIALOG_PASSWORD = "READ_DB_DIALOG_PASSWORD";
  public static String READ_DB_DIALOG_PIPELINE_DB = "READ_DB_DIALOG_PIPELINE_DB";
  public static String READ_DB_DIALOG_PIPELINE_DB_LIST = "READ_DB_DIALOG_PIPELINE_DB_LIST";
  
  public static String LAYOUT_DIALOG = "LAYOUT_DIALOG";
  public static String LAYOUT_DIALOG_ITERATES = "LAYOUT_DIALOG_ITERATES";
  public static String LAYOUT_DIALOG_SPRING_NATURAL_LENGTH = "LAYOUT_DIALOG_SPRING_NATURAL_LENGTH";
  public static String LAYOUT_DIALOG_REPULSION_MULTIPLIER = "LAYOUT_DIALOG_REPULSION_MULTIPLIER";
  public static String LAYOUT_DIALOG_VERTICAL_SPACING = "LAYOUT_DIALOG_VERTICAL_SPACING";
  public static String LAYOUT_DIALOG_HORIZONTAL_SPACING = "LAYOUT_DIALOG_HORIZONTAL_SPACING";
  public static String LAYOUT_DIALOG_FIX_ROOTS = "LAYOUT_DIALOG_FIX_ROOTS";
  public static String LAYOUT_DIALOG_MOVEMENT_LIMIT = "LAYOUT_DIALOG_MOVEMENT_LIMIT";
  public static String LAYOUT_DIALOG_GRAVITY = "LAYOUT_DIALOG_GRAVITY";
  public static String LAYOUT_DIALOG_SHOW_JOB_DETAIL = "LAYOUT_DIALOG_SHOW_JOB_DETAIL";

  public static String CREATE_NODE_DIALOG = "CREATE_NODE_DIALOG";
  public static String CREATE_NODE_DIALOG_NAME = "CREATE_NODE_DIALOG_NAME";
  public static String CREATE_NODE_DIALOG_INPUT_ID_TYPE = "CREATE_NODE_DIALOG_INPUT_ID_TYPE";
  public static String CREATE_NODE_DIALOG_INPUT_ID_TYPE_LIST = "CREATE_NODE_DIALOG_INPUT_ID_TYPE_LIST";
  public static String CREATE_NODE_DIALOG_DB = "CREATE_NODE_DIALOG_DB";
  public static String CREATE_NODE_DIALOG_DB_FILE = "CREATE_NODE_DIALOG_DB_FILE";
  public static String CREATE_NODE_DIALOG_PROGRAM = "CREATE_NODE_DIALOG_PROGRAM";
  public static String CREATE_NODE_DIALOG_PROGRAM_VERSION = "CREATE_NODE_DIALOG_PROGRAM_VERSION";
  public static String CREATE_NODE_DIALOG_PROGRAM_FILE = "CREATE_NODE_DIALOG_PROGRAM_FILE";
  public static String CREATE_NODE_DIALOG_PARAMETERS = "CREATE_NODE_DIALOG_PARAMETERS";
  public static String CREATE_NODE_DIALOG_MODULE = "CREATE_NODE_DIALOG_MODULE";
  public static String CREATE_NODE_DIALOG_MODULE_VERSION = "CREATE_NODE_DIALOG_MODULE_VERSION";
  public static String CREATE_NODE_DIALOG_GFF_SOURCE = "CREATE_NODE_DIALOG_GFF_SOURCE";
  public static String CREATE_NODE_DIALOG_GFF_FEATURE = "CREATE_NODE_DIALOG_GFF_FEATURE";
  
  ModelElement rootElement = new ModelElement(ROOT); 
  
  public PathStepsModel(){
  }

  /**
   * When the model is first created it should contain the Root, and the PathStepsPanel.
  **/
  public void initialise(){
    setRootElement(new ModelElement(ROOT));
    
    //
    //Create a child to contain the PathStepsPanel model
    getRootElement().addChildElement(PATH_STEPS_PANEL, new ModelElement(PATH_STEPS_PANEL));
  }
  
  public ModelElement getRootElement(){
    return rootElement;
  }
  
  public void setRootElement(ModelElement element){
    rootElement = element;
  }
  
  public String toString(){
    return getRootElement().toString();
  }
  
  public void clearMessage() {
    getRootElement().removeProperty(MESSAGE);
  }
  
  public void createMessage(String message) {
    getRootElement().addProperty(MESSAGE, message);
  }
  
}