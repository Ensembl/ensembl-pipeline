package pathsteps.gui;
import java.awt.event.*;

import pathsteps.common.*;

/**
 * Generic class to forward action events onto central handler.
**/
public class KeyEventRouter extends EventRouter implements KeyListener
{
  public KeyEventRouter(Application handler, String key){
    super(handler, key);
  }

  public void keyPressed(KeyEvent e) {
    if(getHandler().getView().getLogger().isLoggingLow()){
      getHandler().getView().getLogger().logLow("KeyEventRouter processing keyPressed(KeyEvent) - passing event with key:"+getKey());
    }
    
    handler.notifyEventForKey(getKey());
  }
  
  public void keyReleased(KeyEvent e) {
  }
  
  public void keyTyped(KeyEvent e) {
  }
  
}
