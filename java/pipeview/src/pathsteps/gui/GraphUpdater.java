package pathsteps.gui;
import javax.swing.*;
import java.awt.event.*;
import org.jgraph.graph.*;
import java.util.*;
import java.awt.*;

/**
 * I use the layout algorithm in LayoutCalculator to change the 
 * input graph model by the input number of iterates. I perform
 * a number of iterations (specified by showInterval)
 * ever time the actionPerformed method is run.
**/
public class GraphUpdater implements ActionListener{
  GraphModel model;
  HashMap graphAttributesByKey;
  HashMap graphCellsByKey;
  int numberOfIterates;
  int showInterval;
  int iteratesPerformed;
  LayoutCalculator calculator;
  int movementLimit;
  int gravity;
  int repulsionMultiplier;
  javax.swing.Timer timer;
  
  public GraphUpdater(
    LayoutCalculator calculator,
    GraphModel model,
    HashMap graphCellsByKey,
    HashMap graphAttributesByKey,
    int numberOfIterates, 
    int showInterval,
    int movementLimit,
    int gravity,
    int repulsionMultiplier
  ){
    this.calculator = calculator;
    this.model = model;
    this.numberOfIterates = numberOfIterates;
    this.showInterval = showInterval;
    this.iteratesPerformed = iteratesPerformed;
    this.graphCellsByKey = graphCellsByKey;
    this.graphAttributesByKey = graphAttributesByKey;
    this.movementLimit = movementLimit;
    this.gravity = gravity;
    this.repulsionMultiplier = repulsionMultiplier;
  }
  
  public void actionPerformed(ActionEvent event){
    //
    //Now read out the Nodes[] array - which will contain the right positions for the nodes.
    int x;
    int y;
    String key;
    ArrayList list = new ArrayList();
    Node[] nodes = calculator.getNodes();
    Map cellAttributeMap;
    Rectangle bounds;
    HashMap cellsToCellAttributeMap = new HashMap();
    DefaultGraphCell cell;

    calculator.layout(
      showInterval, 
      movementLimit, 
      gravity, 
      repulsionMultiplier
    );

    for(int i = 0; i<nodes.length; i++){
      key = nodes[i].lbl;
      x = (new Double(nodes[i].x)).intValue();
      y = (new Double(nodes[i].y)).intValue();
      
      cell = (DefaultGraphCell)getGraphCellsByKey().get(key);
      cellAttributeMap = (Map)getGraphAttributesByKey().get(key);
      bounds = GraphConstants.getBounds(cellAttributeMap);
      bounds.x = x;
      bounds.y = y;
      GraphConstants.setBounds(cellAttributeMap, bounds);
      cellsToCellAttributeMap.put(cell, cellAttributeMap);
    }
    
    getModel().edit(cellsToCellAttributeMap, null, null, null);
    
    numberOfIterates -= showInterval;
    if(numberOfIterates <0){
      timer.stop();

      System.out.println("Afterwards: ---------------");
      for(int i=0; i<nodes.length; i++){
        System.out.println(nodes[i].lbl+" "+nodes[i].x+"-"+nodes[i].y);
      }
      System.out.println("---------------");
    }
  }
  
  private GraphModel getModel(){
    return model;
  }
  
  private int getNumberOfIterates(){
    return numberOfIterates;
  }
  
  private int getShowInterval(){
    return showInterval;
  }
  
  private int getIteratesPerformed(){
    return iteratesPerformed;
  }
  
  private void setIteratesPerformed(int newValue){
    iteratesPerformed = newValue;
  }

  private int getMovementLimit(){
    return movementLimit;
  }
  
  private int getGravity(){
    return gravity;
  }
  
  private int getRepulsionMultiplier(){
    return repulsionMultiplier;
  }
  
  private HashMap getGraphAttributesByKey(){
    return graphAttributesByKey;
  }
  
  private HashMap getGraphCellsByKey(){
    return graphCellsByKey;
  }
}
