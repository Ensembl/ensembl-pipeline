package pathsteps.control;
import pathsteps.common.*;
import pathsteps.model.*;
import pathsteps.gui.*;
import java.util.*;

/**
 * Writes the properties file stored in the layout model elements 
 * PATH_STEPS_PANEL_GRAPH_LAYOUT_CONFIGURATION property onto
 * a file.
**/
public class SaveGraphLayoutConfigurationAction extends AAction{
  public SaveGraphLayoutConfigurationAction(AEventHandler eventHandler){
    super(eventHandler);
  }
  
  protected void doAction(AView view, AModel aModel){
    
    PathStepsModel model = (PathStepsModel)aModel;
    ModelElement panelModel = model.getRootElement().getChildElement(model.PATH_STEPS_PANEL);
    Properties graphLayout = null;
    
    view.update(model);

    
    if(panelModel != null){
      graphLayout = 
        (Properties)panelModel
          .getProperty(
            model.PATH_STEPS_PANEL_GRAPH_LAYOUT_CONFIGURATION
          );
    }else{
      throw new FatalAException("Expecting a model element with key: "+model.PATH_STEPS_PANEL+" - none found");
    }

    if(graphLayout != null){
      view.getApplication().writeGraphLayoutConfiguration(graphLayout);
    }else{
      throw new FatalAException(
        "Expecting a Properties object as a property (keyed by) "+
        model.PATH_STEPS_PANEL_GRAPH_LAYOUT_CONFIGURATION +
        " for the panel model element - none found "
      );
    }
  }
}
