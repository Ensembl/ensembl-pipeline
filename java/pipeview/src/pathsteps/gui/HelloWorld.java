
package pathsteps.gui;

import org.jgraph.JGraph;
import org.jgraph.graph.*;
import javax.swing.BorderFactory;
import java.util.Hashtable;
import java.awt.Rectangle;
import java.awt.Color;
import java.util.Map;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;

public class HelloWorld {

  public ActionListener colorListener;
  public JFrame frame;
  public JScrollPane pane;
  public GraphModel model;
  public DefaultGraphCell cell;
  public DefaultGraphCell cell2;
  public boolean switchMe = false;
  
  public class ColorListener implements ActionListener{
    public void actionPerformed(ActionEvent event){
        Map map = HelloWorld.this.worldAttrib;
        DefaultGraphCell cell = HelloWorld.this.cell;
        java.util.HashMap temp = new java.util.HashMap();
        
        temp.put(cell, map);
        
        temp.put(HelloWorld.this.cell2, HelloWorld.this.helloAttrib);
        
        Rectangle rect = GraphConstants.getBounds(map);
        rect.x +=10;
        rect.y +=10;
        GraphConstants.setBounds(map, rect);
        GraphConstants.setBackground(map, Color.red);
        GraphConstants.setValue(map, "CHANGED");
        
        Rectangle rect2 = GraphConstants.getBounds(HelloWorld.this.helloAttrib);
        rect2.x +=10;
        
        HelloWorld.this.model.edit(temp, null, null, null);
      /*
      if(HelloWorld.this.switchMe){
        System.out.println("Creating a graph based on model :");
        GraphModel model = HelloWorld.this.model;
        JScrollPane pane2 = 
          new JScrollPane(
            new JGraph(HelloWorld.this.model)
          );
        HelloWorld.this.pane = pane2;
        HelloWorld.this.frame.getContentPane().add(pane2);
        HelloWorld.this.frame.invalidate();
        HelloWorld.this.frame.validate();
        HelloWorld.this.frame.repaint();
        HelloWorld.this.switchMe = false;
      }else{
        Map map = HelloWorld.this.worldAttrib;
        DefaultGraphCell cell = HelloWorld.this.cell;
        java.util.HashMap temp = new java.util.HashMap();
        temp.put(cell, map);
        GraphConstants.setBackground(map, Color.red);
        GraphConstants.setValue(map, "CHANGED");
        HelloWorld.this.model.edit(temp, null, null, null);
        HelloWorld.this.frame.getContentPane().remove(HelloWorld.this.pane);
        HelloWorld.this.frame.invalidate();
        HelloWorld.this.frame.validate();
        HelloWorld.this.frame.repaint();
        HelloWorld.this.switchMe = true;
      }
       */
    }
  }

  public HelloWorld(){
    colorListener = new ColorListener();
  }
  
  Map worldAttrib;
  Map helloAttrib;
  
  public static void main(String[] args) {
    
    HelloWorld state = new HelloWorld();

    // Construct Model and Graph
    //
    GraphModel model = new DefaultGraphModel();
    JGraph graph = new JGraph(model);
    //graph.setSelectNewCells(true);

    // Create Nested Map (from Cells to Attributes)
    //
    Map attributes = new Hashtable();

    // Create Hello Vertex
    //
    DefaultGraphCell hello = new DefaultGraphCell("Hello");
    state.cell2 = hello;
    // Create Hello Vertex Attributes
    //
    Map helloAttrib = GraphConstants.createMap();
    attributes.put(hello, helloAttrib);
    state.helloAttrib = helloAttrib;
    
    // Set bounds
    Rectangle helloBounds = new Rectangle(20, 20, 40, 20);
    GraphConstants.setBounds(helloAttrib, helloBounds);
    // Set black border
    GraphConstants.setBorderColor(helloAttrib, Color.black);

    // Add a Port
    //
    DefaultPort hp = new DefaultPort();
    hello.add(hp);

    // Create World Vertex
    //
    DefaultGraphCell world = new DefaultGraphCell("World");
    state.cell = world;

    // Create World Vertex Attributes
    //
    state.worldAttrib = GraphConstants.createMap();
    attributes.put(world, state.worldAttrib);
    // Set bounds
    Rectangle worldBounds= new Rectangle(140, 140, 40, 20);
    GraphConstants.setBounds(state.worldAttrib, worldBounds);
    // Set fill color
    GraphConstants.setBackground(state.worldAttrib, Color.orange);
    GraphConstants.setOpaque(state.worldAttrib, true);
    // Set raised border
    GraphConstants.setBorder(state.worldAttrib, 
           BorderFactory.createRaisedBevelBorder());

    // Add a Port
    //
    DefaultPort wp = new DefaultPort();
    world.add(wp);

    // Create Edge
    //
    DefaultEdge edge = new DefaultEdge();
    
    // Create Edge Attributes
    //
    Map edgeAttrib = GraphConstants.createMap();
    attributes.put(edge, edgeAttrib);
    // Set Arrow
    int arrow = GraphConstants.ARROW_CLASSIC;
    GraphConstants.setLineEnd(edgeAttrib , arrow);
    GraphConstants.setEndFill(edgeAttrib, true);

    // Connect Edge
    //
    ConnectionSet cs = new ConnectionSet(edge, hp, wp);
    Object[] cells = new Object[]{edge, hello, world};

    // Insert into Model
    //
    model.insert(cells, attributes, cs, null, null);

    // Show in Frame
    //
    JFrame frame = new JFrame();
    JButton button = new JButton("Press me");
    button.addActionListener(state.colorListener);

    state.pane = new JScrollPane(graph);
    frame.getContentPane().add(state.pane, BorderLayout.CENTER);
    frame.getContentPane().add(button, BorderLayout.SOUTH);
    frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    frame.pack();
    frame.setVisible(true);
    state.frame = frame;
    state.model = model;
  }  
}
