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
    if(panelModel != null){
      graphLayout = 
        (Properties)panelModel
          .getChildElement(model.PATH_STEPS_PANEL)
          .getProperty(model.PATH_STEPS_PANEL_GRAPH_LAYOUT_CONFIGURATION);
    }

    if(graphLayout != null){
      view.getApplication().writeGraphLayoutConfiguration(graphLayout);
    }
  }
}
