package pathsteps.gui;
import java.awt.event.*;

import pathsteps.common.*;

/**
 * Generic class to forward action events onto central handler.
**/
public class WindowCloseEventRouter extends EventRouter implements WindowListener
{
  public WindowCloseEventRouter(Application handler, String key){
    super(handler, key);
  }
  
  public void windowActivated(java.awt.event.WindowEvent windowEvent) {
  }
  
  public void windowClosed(java.awt.event.WindowEvent windowEvent) {
  }
  
  public void windowClosing(java.awt.event.WindowEvent windowEvent) {
    if(getHandler().getView().getLogger().isLoggingLow()){
      getHandler().getView().getLogger().logLow("WindowCloseEventRouter processing windowClosing(WindowEvent) - passing event with key:"+getKey());
    }
    getHandler().notifyEventForKey(getKey());
  }
  
  public void windowDeactivated(java.awt.event.WindowEvent windowEvent) {
  }
  
  public void windowDeiconified(java.awt.event.WindowEvent windowEvent) {
  }
  
  public void windowIconified(java.awt.event.WindowEvent windowEvent) {
  }
  
  public void windowOpened(java.awt.event.WindowEvent windowEvent) {
  }
}
