package pathsteps.gui;

import java.awt.event.*;
import javax.swing.*;

import pathsteps.common.*;

/**
 * Generic class to forward mouse left click events onto central handler.
**/
public class 
  MouseLeftClickEventRouter 
extends 
  EventRouter
implements 
  MouseListener
{
 
  public MouseLeftClickEventRouter(Application handler, String key){
    super(handler, key);
  }
  
  public void mouseClicked(MouseEvent mouseEvent) {
    if(SwingUtilities.isLeftMouseButton(mouseEvent)){
      getHandler().notifyEventForKey(key);
    }
  }
  
  public void mouseEntered(MouseEvent mouseEvent) {}
  
  public void mouseExited(MouseEvent mouseEvent) {}
  
  public void mousePressed(MouseEvent mouseEvent) {}
  
  public void mouseReleased(MouseEvent mouseEvent) {}
}
