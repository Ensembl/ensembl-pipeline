package pathsteps.test;
import pathsteps.model.*;

import java.util.*;
import junit.framework.*;

public class CreateNewLayoutConfigurationDialogTest extends TestCase{
  TestRunner runner;
  
  public CreateNewLayoutConfigurationDialogTest(TestRunner runner){
    super("testAction");
    this.runner = runner;
  }
  
  public void setUp(){}
  
  public void tearDown(){}
  
  public void testAction(){
    
    String iterates = getStringFromHistoryOrConfig(getRunner(), PathStepsModel.LAYOUT_DIALOG_ITERATES);
    String gravity = getStringFromHistoryOrConfig(getRunner(), PathStepsModel.LAYOUT_DIALOG_GRAVITY);
    String horizontal_space = getStringFromHistoryOrConfig(getRunner(), PathStepsModel.LAYOUT_DIALOG_HORIZONTAL_SPACING);
    String vertical_space = getStringFromHistoryOrConfig(getRunner(), PathStepsModel.LAYOUT_DIALOG_VERTICAL_SPACING);
    String limit = getStringFromHistoryOrConfig(getRunner(), PathStepsModel.LAYOUT_DIALOG_MOVEMENT_LIMIT);
    String repulsion = getStringFromHistoryOrConfig(getRunner(), PathStepsModel.LAYOUT_DIALOG_REPULSION_MULTIPLIER);
    String spring_length = getStringFromHistoryOrConfig(getRunner(), PathStepsModel.LAYOUT_DIALOG_SPRING_NATURAL_LENGTH);
    
    PathStepsModel model = getRunner().getModel();
    
    getRunner().getApplication().notifyEventForKey(getRunner().getView().SET_NEW_LAYOUT_PREFERENCES_KEY);

    Properties layoutDialogProperties = getRunner().getView().getLayoutDialogData();
    
    assertEquals(iterates, layoutDialogProperties.getProperty(PathStepsModel.LAYOUT_DIALOG_ITERATES));
    assertEquals(gravity, layoutDialogProperties.getProperty(PathStepsModel.LAYOUT_DIALOG_GRAVITY));
    assertEquals(horizontal_space, layoutDialogProperties.getProperty(PathStepsModel.LAYOUT_DIALOG_HORIZONTAL_SPACING));
    assertEquals(vertical_space, layoutDialogProperties.getProperty(PathStepsModel.LAYOUT_DIALOG_VERTICAL_SPACING));
    assertEquals(limit, layoutDialogProperties.getProperty(PathStepsModel.LAYOUT_DIALOG_MOVEMENT_LIMIT));
    assertEquals(repulsion, layoutDialogProperties.getProperty(PathStepsModel.LAYOUT_DIALOG_REPULSION_MULTIPLIER));
    assertEquals(spring_length, layoutDialogProperties.getProperty(PathStepsModel.LAYOUT_DIALOG_SPRING_NATURAL_LENGTH));
  }
  
  private String getStringFromHistoryOrConfig(TestRunner runner, String key){

    Properties configuration = runner.getApplication().getConfiguration();
    Properties history = runner.getApplication().getHistory();
    String returnString = history.getProperty(key);
    
    if(returnString == null || returnString.trim().length() <=0){
      returnString = configuration.getProperty(key);
    }
    
    if(returnString == null || returnString.trim().length() <=0){
      fail("Could not get value for: "+key+" from history or contig");
    }

    return returnString;
  }
  
  public TestRunner getRunner(){
    return runner;
  }
}
