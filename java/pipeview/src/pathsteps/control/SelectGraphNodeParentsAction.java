package pathsteps.control;
import pathsteps.common.*;
import pathsteps.model.*;

import java.util.*;

public class SelectGraphNodeParentsAction extends AAction{
  
  public SelectGraphNodeParentsAction(AEventHandler eventHandler){
    super(eventHandler);
  }

  protected void doAction(AView view, AModel aModel){
    PathStepsModel model = (PathStepsModel)aModel;
    ModelElement rootElement;
    ModelElement selectedElement;
    pathsteps.gui.SetUtil utils = new pathsteps.gui.SetUtil(getLogger());
    String nodeName = view.getNameOfSelectedNode();

    if(nodeName != null){
      rootElement = model.getRootElement().getChildElement(model.PATH_STEPS_PANEL);
      selectedElement = utils.findRelatedElementWithName(rootElement, nodeName);
      utils.clearPropertyFromChildrenOfNode(rootElement, ModelElement.SELECTED);
      utils.addPropertyToParentsOfNode(selectedElement, rootElement, ModelElement.SELECTED, Boolean.TRUE.toString());
    }
  }
}
