package pathsteps.control;
import pathsteps.common.*;

public class RootFrameCloseAction extends AAction{
  public RootFrameCloseAction(AEventHandler eventHandler){
    super(eventHandler);
  }
  
  protected void doAction(AView view, AModel model){
    view.shutDown();
  }
}
