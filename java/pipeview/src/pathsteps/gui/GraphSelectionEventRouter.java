package pathsteps.gui;
import java.awt.event.*;
import org.jgraph.event.*;

import pathsteps.common.*;

/**
 * Generic class to forward GraphSelectionEvents onto central handler.
**/
public class GraphSelectionEventRouter extends EventRouter implements GraphSelectionListener
{
  public GraphSelectionEventRouter(Application handler, String key){
    super(handler, key);
  }
  
  public void valueChanged(org.jgraph.event.GraphSelectionEvent graphSelectionEvent) {
    if(getHandler().getView().getLogger().isLoggingLow()){
      getHandler().getView().getLogger().logLow("GraphSelectionEventRouter processing valueChanged(GraphSelectionEvent) - passing event with key:"+getKey());
    }
    getHandler().notifyEventForKey(getKey());    
  }
}
