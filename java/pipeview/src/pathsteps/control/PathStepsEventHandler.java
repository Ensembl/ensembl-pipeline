package pathsteps.control;
import java.util.*;
import javax.swing.*;
import java.awt.event.*;
import pathsteps.common.*;

public class PathStepsEventHandler implements AEventHandler{
  HashMap actions;
  private Application application;
  private ALogger logger;
  
  public PathStepsEventHandler(){
    setActions(new HashMap());
  }
  
  public void initialise(){
    if(getApplication() == null){
      throw new FatalAException("Attempt to initialise PathStepsEventHandler without first defining Application");
    }

    getActions().put(AView.SET_NEW_LAYOUT_PREFERENCES_KEY, new CreateNewLayoutConfigurationDialogAction(this));
    getActions().put(AView.READ_NEW_PIPELINE_DB_KEY, new CreateNewReadPipelineDBDialogAction(this));
    getActions().put(AView.CLOSE_ROOT_FRAME_KEY, new RootFrameCloseAction(this));
    getActions().put(AView.DATABASE_DATA_CHANGE_KEY, new DatabaseDataChangeAction(this));
    getActions().put(AView.FIND_ALL_DATABASES_KEY, new FindAllDatabasesAction(this));
    getActions().put(AView.CANCEL_READ_KEY, new CancelReadAction(this));
    getActions().put(AView.READ_SPECIFIED_DB_KEY, new ReadSpecifiedDBAction(this));
    getActions().put(AView.REFRESH_BUTTON_KEY, new RereadSpecifiedDBAction(this));
    getActions().put(AView.LAYOUT_BUTTON_KEY, new ApplyGraphLayoutAction(this));
    getActions().put(AView.SET_LAYOUT_PREFERENCES_KEY, new SetLayoutConfigurationAction(this));
    getActions().put(AView.CANCEL_LAYOUT_PREFERENCES_KEY, new CancelLayoutConfigurationAction(this));
    getActions().put(AView.SELECT_GRAPH_NODE_PARENTS_KEY, new SelectGraphNodeParentsAction(this));
    getActions().put(AView.CREATE_NODE_DIALOG_KEY, new CreateNodeDialogAction(this));
    getActions().put(AView.SAVE_GRAPH_LAYOUT_CONFIGURATION_KEY, new SaveGraphLayoutConfigurationAction(this));
  }

  
  public void notifyEventForKey(String key){
    AModel model = getApplication().getModel();
    AView view = getApplication().getView();
    ALogger eventLogger = getLogger();
    
    AAction action = (AAction)getActions().get(key);
    
    if(action == null){
      throw new FatalAException("No action found for key: "+key);
    }
    
    try{
      
      if(eventLogger.isLoggingLow()){
        eventLogger.logLow("Pre-action update for action: "+key);
      }
      
      model.clearMessage();
      view.showMessage("");
      
      view.update(model);
      
      if(eventLogger.isLoggingLow()){
        eventLogger.logLow("Running action with key: "+key);
      }
      
      action.processAction();
      
      if(eventLogger.isLoggingLow()){
        eventLogger.logLow("After notifying EventHandler for action: "+key);
      }
      
      view.read(model);
      
      if(eventLogger.isLoggingLow()){
        eventLogger.logLow("After post-action read for action: "+key);
      }
      
    }catch(NonFatalAException exception){
      
      view.showMessage(exception.getMessage());
      model.createMessage(exception.getMessage());
      view.read(model);
    
      if(eventLogger.isLoggingLow()){
        eventLogger.logLow("NonFatalProblem ProcessAction: exception"+exception.getMessage(), exception);
        if(eventLogger.isLoggingMedium()){
          eventLogger.logMedium("ProcessAction: exception"+exception.getMessage(), exception);
        }
      }
      
    }catch(Throwable exception){
      
      String message="";
      if(exception.getMessage() != null){
        message = exception.getMessage();
      }
      
      if(eventLogger.isLoggingLow()){
        eventLogger.logLow("General throwable: exception"+message, exception);
        eventLogger.logLow("Showing fatal dialog");
        eventLogger.logLow("Shutting application down");
      }
      
      view.showFatalDialog(exception.getClass()+" "+message);
      view.shutDown();
    }

    if(eventLogger.isLoggingLow()){
      eventLogger.logLow("Finished Notifying EventHandler for action with key: "+key);
    }
  }
  
  public void doActionForKey(String key){
    notifyEventForKey(key);
  }
  
  private HashMap getActions(){
    return actions;
  }
  
  private void setActions(HashMap actions){
    this.actions = actions;
  }
  
  public void loadActions(){
    //getActions().add("FAKE_ACTION_KEY", new FakeAction());
  }
  
  public ALogger getLogger(){
    return logger;
  }
  
  public void setLogger(ALogger logger){
    this.logger = logger;
  }
  
  public Application getApplication(){
    return application;
  }
  
  public void setApplication(Application application){
    this.application = application;
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
