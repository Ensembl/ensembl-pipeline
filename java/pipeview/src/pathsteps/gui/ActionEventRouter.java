package pathsteps.gui;
import java.awt.event.*;

import pathsteps.common.*;

/**
 * Generic class to forward action events onto central handler.
**/
public class ActionEventRouter 
extends EventRouter
implements ActionListener
{
  public ActionEventRouter(Application handler, String key){
    super(handler, key);
  }
  
  public void actionPerformed(ActionEvent event){
    if(getHandler().getView().getLogger().isLoggingLow()){
      getHandler().getView().getLogger().logLow("ActionEventRouter processing actionPerformed(ActionEvent) - passing event with key:"+getKey());
    }
    getHandler().notifyEventForKey(getKey());
  }
}
