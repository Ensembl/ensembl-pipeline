package pathsteps.gui;
import java.awt.event.*;

import pathsteps.common.*;

/**
 * Generic class to forward any old events onto central handler.
**/
public abstract class EventRouter
{
  String key;
  Application handler;
  
  public EventRouter(Application handler, String key){
    this.handler = handler;
    this.key = key;
  }
  
  protected Application getHandler(){
    return handler;
  }
  
  protected String getKey(){
    return key;
  }
}
