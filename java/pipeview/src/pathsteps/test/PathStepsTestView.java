package pathsteps.test;
import pathsteps.common.*;
import pathsteps.model.*;

import java.util.*;

/**
 * Swing implementation of the AView interface
**/
public class PathStepsTestView implements AView
{
  private Application application;
  private String messageText;
  private PathStepsModel model;
  
  private boolean layoutConfigurationDialogOpen = false;
  
  public Application getApplication(){
    return application;
  }
  
  public void setApplication(Application application){
    this.application = application;
  }
  
  public static PathStepsTestView create(){
    return new PathStepsTestView();
  }
  
  public void read(AModel model){
    this.model = (PathStepsModel)model;
    System.out.println("Test view has done a read - new model is:");
    System.out.println(model.toString());
  }
  
  public PathStepsModel getModel(){
    return model;
  }
  
  public void update(AModel model) {
  }
  
  public void initialise(){
  }
  
  public void show(){
  }
  
  public void openReadPipelineDBDialog(){
    
  }

  public void showMessage(String message){
    messageText = message;
  }
  

  public void showFatalDialog(String message){
    messageText = message;
  }
  
  public void shutDown(){
  }
  
  public void closeReadPipelineDBDialog() {
  }
  
  public ALogger getLogger() {
    return null;
  }
  
  public void setLogger(ALogger logger) {
  }
  
  public boolean isReadPipelineDBDialogOpen() {
    return false;
  }
  
  public void requestFocusOnPipelineDBDialog() {
  }

  public void openLayoutConfigurationDialog(){
  }

  public void closeLayoutConfigurationDialog(){
  }
  
  public boolean isLayoutConfigurationDialogOpen() {
    return layoutConfigurationDialogOpen;
  }
  
  public void requestFocusOnLayoutConfigurationDialog() {
  }
  
  public void closeCreateNodeDialog() {
  }
  
  public boolean isCreateNodeDialogOpen() {
    return false;
  }
  
  public void openCreateNodeDialog() {
    
  }
  
  public void requestFocusOnCreateNodeDialog() {
  }
  
  public void reshow() {
  }
  
  public String getNameOfSelectedNode() {
    return null;
  }

  /**
   * Specific to this TEST view -to help internal testing.
  **/
  public Properties getLayoutDialogData(){
    Properties dialogProps = new Properties();
    ModelElement element = getModel().getRootElement().getChildElement(PathStepsModel.LAYOUT_DIALOG);
    element = getModel().getRootElement().getChildElement(PathStepsModel.LAYOUT_DIALOG);
    if(element == null){
      throw new FatalAException("Should have found a non-null element in RootElement with key: "+PathStepsModel.LAYOUT_DIALOG);
    }

    dialogProps.put(
      PathStepsModel.LAYOUT_DIALOG_ITERATES, 
      element.getProperty(PathStepsModel.LAYOUT_DIALOG_ITERATES)
    );
    
    dialogProps.put(
      PathStepsModel.LAYOUT_DIALOG_GRAVITY, 
      element.getProperty(PathStepsModel.LAYOUT_DIALOG_GRAVITY)
    );
    dialogProps.put(
      PathStepsModel.LAYOUT_DIALOG_HORIZONTAL_SPACING, 
      element.getProperty(PathStepsModel.LAYOUT_DIALOG_HORIZONTAL_SPACING)
    );
    dialogProps.put(
      PathStepsModel.LAYOUT_DIALOG_VERTICAL_SPACING, 
      element.getProperty(PathStepsModel.LAYOUT_DIALOG_VERTICAL_SPACING)
    );
    dialogProps.put(
      PathStepsModel.LAYOUT_DIALOG_MOVEMENT_LIMIT, 
      element.getProperty(PathStepsModel.LAYOUT_DIALOG_MOVEMENT_LIMIT)
    );
    dialogProps.put(
      PathStepsModel.LAYOUT_DIALOG_REPULSION_MULTIPLIER, 
      element.getProperty(PathStepsModel.LAYOUT_DIALOG_REPULSION_MULTIPLIER)
    );
    dialogProps.put(
      PathStepsModel.LAYOUT_DIALOG_SPRING_NATURAL_LENGTH, 
      element.getProperty(PathStepsModel.LAYOUT_DIALOG_SPRING_NATURAL_LENGTH)
    );
    
    return dialogProps;
  }
  
  public void setLayoutConfigurationDialogOpen(boolean value) {
    layoutConfigurationDialogOpen = value;
  }
  
  public void applyGraphLayout() {
  }

  public boolean isLoggingLow(){
    return true;
  }
  
  public void logLow(String message){
    
  }
  
  public void logLow(String message, Throwable exception) {
    
  }
  
  public boolean isLoggingMedium() {
    return true;
  }
  
  public void logMedium(String message) {
    
  }
  
  public void logMedium(String message, Throwable exception) {
    
  }
  
  public boolean isLoggingHigh(){
    return true;
  }
  
  public void logHigh(String message) {
    
  }
  
  public void logHigh(String message, Throwable exception) {
    
  }  
}