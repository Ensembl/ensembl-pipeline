package pathsteps.gui.graph;

import org.jgraph.*;
import org.jgraph.graph.*;
import java.util.*;
import javax.swing.*;
import javax.swing.tree.*;
import pathsteps.model.*;

/**
 * The point of the this class is to have an overridden createVertexView method,
 * so we can create our own CellView class, (PipeViewCell), which in turn
 * can render more 'stuff' (via the PipeViewCellRenderer) than the default cell 
 * renderer.
**/
public class PipeViewGraph extends JGraph{

  public PipeViewGraph(GraphModel model){
    super(model);
  }
/*  
  protected VertexView createVertexView(
    Object cell,
    GraphModel model,
    CellMapper cellMapper
  ){
    return new PipeViewView(cell, model, cellMapper);
  }
*/
  public String convertValueToString(Object value) {
 		CellView view = (CellView) value;
		Object newValue;
    Object property;
    Map statusMap;
    String finishedString = "";
    String createdString = "";
    String failedString = "";
    String submittedString = "";
    String readingString = "";
    String runningString = "";
    String voidString = "";
    boolean showDetail = false;
    java.util.Map attributes = view.getAllAttributes();
    
    if (view != null) {
      if(view instanceof VertexView){

        property = attributes.get(ModelElement.SHOW_DETAIL);
        showDetail = Boolean.valueOf((String)property).booleanValue();
        if(showDetail){
          attributes = view.getAllAttributes();
          property = attributes.get(ModelElement.FINISHED_COUNT);
          if(property != null){
            finishedString = ((Integer)property).toString();
          }

          property = attributes.get(ModelElement.ANALYSIS_STATUS_COUNT);
          if(property != null){
            statusMap = (Map)property;
            if(statusMap.get(ModelElement.CREATED) != null){
              createdString = ((Integer)statusMap.get(ModelElement.CREATED)).toString(); 
            }

            if(statusMap.get(ModelElement.RUNNING) != null){
              runningString = ((Integer)statusMap.get(ModelElement.RUNNING)).toString(); 
            }

            if(statusMap.get(ModelElement.FAILED) != null){
              failedString = ((Integer)statusMap.get(ModelElement.FAILED)).toString(); 
            }
          }
        }
      }
      
      newValue = GraphConstants.getValue(view.getAllAttributes());
      
      if (newValue != null){
        if(showDetail){
          value = 
            "<html> "+
            "<strong><em>" +newValue+"</em></strong>"+
            " <br>FIN: "+finishedString+
            " <br>CRE: "+createdString+
            " <br>RUN: "+runningString+
            " <br>FAI: "+failedString+
            " <html>";
        }else{
          value = newValue;
        }
      }else{
				value = view.getCell();
      }
		}
    
		if (
      value instanceof DefaultMutableTreeNode
			&& ((DefaultMutableTreeNode) value).getUserObject() != null
    ){
      
			return ((DefaultMutableTreeNode) value).getUserObject().toString();
      
    } else if (value != null){
      
			return value.toString();
      
    }
    
		return null;
	}

}
