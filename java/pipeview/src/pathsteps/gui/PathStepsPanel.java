package pathsteps.gui;
import pathsteps.gui.graph.PipeViewGraph;
import pathsteps.model.*;
import pathsteps.common.*;
import javax.swing.*;
import java.awt.*;
import java.awt.geom.*;
import java.util.*;

import org.jgraph.*;
import org.jgraph.graph.*;
import javax.swing.JFrame;
import javax.swing.JScrollPane;
import javax.swing.BorderFactory;
import java.util.Hashtable;
import java.awt.Rectangle;
import java.awt.Color;
import java.util.Map;

public class PathStepsPanel extends JPanel{
  public static int ROW_SPACING = 100; //pixels
  public static int COLUMN_SPACING = 100; //pixels
  public static int RECTANGLE_WIDTH = 125;
  public static int RECTANGLE_HEIGHT = 75;
  //stores the maps of modelelementKey -> (x-y) position of each analysis name
  private HashMap keyToPositionMap = new HashMap();
  
  PathStepsSwingView view;  
  PathStepsModel model;
  
  private GraphModel graphModel;
  private JGraph graph;
  private HashMap graphAttributeMap = new HashMap();
  private HashMap cellMap = new HashMap();
  private HashMap elementPortMap = new HashMap();
  private HashSet drawnSet = new HashSet();
  private JScrollPane scrollPane;
  private FontMetrics fontMetrics;

  private Node[] nodes;
  Edge[] edges;
  
  public class Link{
    public String from;
    public String to;
  }

  /**
   * Orders elements in reverse order of number of children - those with
   * lots of kids go first
  **/
  public class InverseChildNumberComparator implements Comparator{
    
    public int compare(Object node1, Object node2) {
      ModelElement element1 = (ModelElement)node1;
      ModelElement element2 = (ModelElement)node2;
      
      return element2.getChildElements().size() - element1.getChildElements().size();
    }
    
    public boolean equals(Object node1, Object node2){
      ModelElement element1 = (ModelElement)node1;
      ModelElement element2 = (ModelElement)node2;
      
      return (element1.getChildElements().size() == element2.getChildElements().size());
    }
  }
  
  public PathStepsPanel(PathStepsSwingView view){
    super(true);
    this.view = view;
    setLayout(new GridBagLayout());
    setGraphModel(new DefaultGraphModel());
    setGraph(new PipeViewGraph(getGraphModel()));
    setScrollPane(new JScrollPane(getGraph()));
    getScrollPane().setPreferredSize(new Dimension(700, 700));
    getScrollPane().setMinimumSize(new Dimension(700, 700));
    GridBagConstraints constraints = new GridBagConstraints();
    constraints.gridx = 0;
    constraints.gridy = 0;
    constraints.weightx= 1;
    constraints.weighty = 1;
    constraints.fill = GridBagConstraints.BOTH;
    add(getScrollPane(), constraints);
  }
  
  private int findNode(String lable){
    for (int i = 0 ; i < nodes.length ; i++) {
      if (nodes[i].lbl.equals(lable)) {
        return i;
      }
    }
    return -10;
  }
  
  public void read(PathStepsModel model){
    int width = 0;
    int height = 0;
    String temp = null;
    int maxWidth = 0;
    int maxHeight = 0;
    int minWidth = 0;
    int minHeight = 0;
    
    if(getView().getLogger().isLoggingMedium()){
      getView().getLogger().logMedium("Started read");
    }
    setModel(model);
    createGraphModel(model);
    recreateGraph();
    
    width = (new Double(getGraph().getPreferredScrollableViewportSize().getWidth())).intValue();
    height = (new Double(getGraph().getPreferredScrollableViewportSize().getHeight())).intValue();
    
    temp = 
      getView().getApplication().getConfiguration().getProperty(
        PathStepsSwingView.MAXIMUM_SCROLLPANE_WIDTH
      );
    
    maxWidth = 0;
    if(temp !=null && temp.trim().length()>0){
      maxWidth = Integer.valueOf(temp).intValue();
    }else{
      maxWidth = 1000;
      if(getView().getLogger().isLoggingHigh()){
        getView().getLogger().logHigh("Could not find config property MAXIMUM_SCROLLPANE_WIDTH - using 1000");
      }
    }

    temp = getView().getApplication().getConfiguration().getProperty(PathStepsSwingView.MAXIMUM_SCROLLPANE_HEIGHT);
    if(temp !=null && temp.trim().length()>0){
      maxHeight = Integer.valueOf(temp).intValue();
    }else{
      maxHeight = 1000;
      if(getView().getLogger().isLoggingHigh()){
        getView().getLogger().logHigh("Could not find config property MAXIMUM_SCROLLPANE_HEIGHT - using 1000");
      }
    }
    
    temp = getView().getApplication().getConfiguration().getProperty(PathStepsSwingView.MINIMUM_SCROLLPANE_HEIGHT);
    if(temp !=null && temp.trim().length()>0){
      minHeight = Integer.valueOf(temp).intValue();
    }else{
      minHeight = 1000;
      if(getView().getLogger().isLoggingHigh()){
        getView().getLogger().logHigh("Could not find config property MINIMUM_SCROLLPANE_HEIGHT - using 1000");
      }
    }

    temp = getView().getApplication().getConfiguration().getProperty(PathStepsSwingView.MINIMUM_SCROLLPANE_WIDTH);
    if(temp !=null && temp.trim().length()>0){
      minWidth = Integer.valueOf(temp).intValue();
    }else{
      minWidth = 1000;
      if(getView().getLogger().isLoggingHigh()){
        getView().getLogger().logHigh("Could not find config property MINIMUM_SCROLLPANE_WIDTH - using 1000");
      }
    }

    if(width < minWidth){
      width = minWidth;
    }
    
    if(width > maxWidth){
      width = maxWidth;
    }

    if(height < minHeight){
      height = minHeight;
    }
    
    if(height > maxHeight){
      height = maxHeight;
    }


    getScrollPane().setMinimumSize(new Dimension(width, height));
    getScrollPane().setPreferredSize(new Dimension(width, height));
    setMinimumSize(new Dimension(width, height));
    setPreferredSize(new Dimension(width, height));
    
    if(getView().getLogger().isLoggingHigh()){
      getView().getLogger().logHigh(
        "The scrollpane's preferred viewport size is: "+
        getGraph().getPreferredScrollableViewportSize().getWidth() +
        ", "+getGraph().getPreferredScrollableViewportSize().getHeight()
      );
      getView().getLogger().logHigh(
        "The preferred min size is set to: "+getMinimumSize().width+", "+getMinimumSize().height
      );
    }
    
    if(getView().getLogger().isLoggingMedium()){
      getView().getLogger().logMedium("Finished read");
    }
  }
  
  private void createGraphModel(PathStepsModel model){
    java.util.List list = new ArrayList();
    ModelElement rootElement = null;
    ModelElement layoutDialogElement = null;
    int row = 0;
    Object[] roots;
    Comparator inverseChildNumber =  new PathStepsPanel.InverseChildNumberComparator();
    
    int xSeparation;
    int ySeparation;
    int naturalEdgeLength;
    int repulsionMultiplier;
    int iterateNumber;
    int movementLimit;
    int gravity;
    boolean fixRoots;
    Properties positionProperties;

    if(getView().getLogger().isLoggingHigh()){
      getView().getLogger().logHigh("Started creating graph");
    }

    if(model != null){
      rootElement = model.getRootElement().getChildElement(PathStepsModel.PATH_STEPS_PANEL);
      if(rootElement != null){
        if(rootElement.isNew()){

          layoutDialogElement = model.getRootElement().getChildElement(PathStepsModel.LAYOUT_DIALOG);
          xSeparation = getNumberFromModel(layoutDialogElement, model.LAYOUT_DIALOG_HORIZONTAL_SPACING);
          ySeparation = getNumberFromModel(layoutDialogElement, model.LAYOUT_DIALOG_VERTICAL_SPACING);
          naturalEdgeLength = getNumberFromModel(layoutDialogElement, model.LAYOUT_DIALOG_SPRING_NATURAL_LENGTH);
          repulsionMultiplier = getNumberFromModel(layoutDialogElement, model.LAYOUT_DIALOG_REPULSION_MULTIPLIER);
          iterateNumber = getNumberFromModel(layoutDialogElement, model.LAYOUT_DIALOG_ITERATES);
          movementLimit = getNumberFromModel(layoutDialogElement, model.LAYOUT_DIALOG_MOVEMENT_LIMIT);
          gravity = getNumberFromModel(layoutDialogElement, model.LAYOUT_DIALOG_GRAVITY);
          fixRoots = getBooleanFromModel(layoutDialogElement, model.LAYOUT_DIALOG_FIX_ROOTS);

          getKeyToPositionMap().clear();
          getGraphAttributeMap().clear();
          getCellMap().clear();
          getKeyToPositionMap().clear();
          getDrawnSet().clear();
          setGraphModel(new DefaultGraphModel());

          list.add(rootElement);
          list = listChildren(list);
          Collections.sort(list, inverseChildNumber);
          //we want to rearrange this list according to which connected components each of
          //these elements belong to.
          SetUtil util = new SetUtil(getView().getLogger());
          list = util.sortRootsByConnectivity(list);

          while(!list.isEmpty()){
            storeKeysAndPositions(list,row);
            list = listChildren(list);
            row++;
          }

          Set keySet = getKeyToPositionMap().keySet();
          Iterator keys = keySet.iterator();
          String key = null;
          Node node = null;
          Dyad position = null;
          list.clear();
          list.add(rootElement);
          list = listChildren(list);
          
          //list contains the first row - the 'roots'. store their keys.
          ArrayList rootKeys = new ArrayList();
          for(int i=0; i<list.size(); i++){
            rootKeys.add(((ModelElement)list.get(i)).getKey());
          }

          int counter=0;
          nodes = new Node[keySet.size()];
          while(keys.hasNext()){
            key = (String)keys.next();
            node = new Node();
            positionProperties = 
              (Properties)
              rootElement.getProperty(model.PATH_STEPS_PANEL_GRAPH_LAYOUT_CONFIGURATION);
            
            if(positionProperties.getProperty(key)!= null){
              position = parsePositionString(positionProperties.getProperty(key));
              node.x = position.x();
              node.y = position.y();
            }else{
              position = (Dyad)getKeyToPositionMap().get(key);
              node.x = xSeparation * position.x() + 10;
              node.y = ySeparation * position.y() + 10;
            }
            
            node.lbl = key;
            if(getView().getLogger().isLoggingHigh()){
              getView().getLogger().logHigh("node: "+key+":"+node.x+"-"+node.y);
            }

            //make fixed if root.
            if(fixRoots && rootKeys.contains(key)){
              if(getView().getLogger().isLoggingHigh()){
                getView().getLogger().logHigh("FIXING node with key: "+key);
              }
              node.fixed = true;
            }
            nodes[counter]=node;
            counter++;
          }

          createGraphModelNodesAndEdges(rootElement, nodes);

          rootElement.age(); //we've read the panel. Don't re-read till a db fetch has been done
          
          recreateGraph();

          
        }else{

          applyEditsFromChangedModelToGraphModel();
          
        }
        
      }else{
        if(getView().getLogger().isLoggingHigh()){
          getView().getLogger().logHigh("Path steps panel model is null - not updating view");
        }
      }
    }else{
      if(getView().getLogger().isLoggingHigh()){
        getView().getLogger().logHigh("Model is null - created no cells");
      }
    }
  }

  private void applyEditsFromChangedModelToGraphModel(){
    ModelElement rootElement = getModel().getRootElement().getChildElement(getModel().PATH_STEPS_PANEL);
    if(getView().getLogger().isLoggingHigh()){
      getView().getLogger().logHigh("UPDATING existing GraphModel model, NOT recreating");
    }

    //
    //make a pass through the model, and set attributes on the 'selected' nodes, so the JGraph will
    //take up the new attributes.
    SetUtil utils = new SetUtil(getView().getLogger());
    Collection allNodes = new HashSet();
    Iterator nodes;
    ModelElement node;
    boolean selected;
    Map cellsToAttributeMap = new HashMap();
    Map cellAttributes;
    DefaultGraphCell cell;

    allNodes = (Collection)rootElement.getProperty(model.PATH_STEPS_PANEL_ALL_NODES);
    nodes = allNodes.iterator();
    while(nodes.hasNext()){
      node = (ModelElement)nodes.next();
      cellAttributes = (Map)getGraphAttributeMap().get(node.getKey());
      cell = (DefaultGraphCell)getCellMap().get(node.getKey());
      if(cellAttributes == null){
        throw new FatalAException("Cannot find the attribute map for cell of key "+node.getKey());
      }

      if(
        node.getProperty(ModelElement.SELECTED) != null &&
        node.getProperty(ModelElement.SELECTED).equals(Boolean.TRUE.toString())
      ){
        if(getView().getLogger().isLoggingMedium()){
          getView().getLogger().logMedium("Setting node "+node.getKey()+" to selected");
        }

        GraphConstants.setBorderColor(cellAttributes, Color.red);
      }else{
        GraphConstants.setBorderColor(cellAttributes, Color.black);
      }
      cellsToAttributeMap.put(cell, cellAttributes);
    }

    getGraphModel().edit(cellsToAttributeMap, null, null, null);
  }
  
  private void createGraphModelNodesAndEdges(ModelElement graphRoot, Node[] nodePositionsAndLabels){
    //
    //Now read out the Nodes[] array - which will contain the right positions for the nodes.
    int width;
    int row;
    ArrayList list = new ArrayList();
    Comparator inverseChildNumber = new PathStepsPanel.InverseChildNumberComparator();
    Map elementMap = (Map)graphRoot.getProperty(PathStepsModel.PATH_STEPS_PANEL_ALL_NODES_MAP);
    String showDetail = (String)graphRoot.getProperty(PathStepsModel.PATH_STEPS_PANEL_SHOW_JOB_DETAIL);
    String finalLabel = "";
    String minimumString = "SubmitSlice";
    int minimumWidth = getFontMetrics().stringWidth(minimumString)+15;
    int rectangleHeight = RECTANGLE_HEIGHT;
    
    if(!Boolean.valueOf(showDetail).booleanValue()){
      rectangleHeight = RECTANGLE_HEIGHT/4;
    }
    
    for(int i=0; i<nodes.length; i++){
      
      width = getFontMetrics().stringWidth(nodePositionsAndLabels[i].lbl)+15;
      if(width < minimumWidth){
        width = minimumWidth;
      }
      
      drawRectangle(
        (ModelElement)elementMap.get(nodePositionsAndLabels[i].lbl),
        (new Long(Math.round(nodes[i].x))).intValue(), 
        (new Long(Math.round(nodes[i].y))).intValue(), 
        width, 
        rectangleHeight, 
        nodePositionsAndLabels[i].lbl
      );
    }

    row=0;
    list.clear();
    list.add(graphRoot);
    list = listChildren(list);
    Collections.sort(list, inverseChildNumber);
    while(!list.isEmpty()){
      drawLinks(row, list);
      Collections.sort(list, inverseChildNumber);
      list = listChildren(list);
      row++;
    }
  }
  
  private void recreateGraph(){
    remove(getScrollPane());
    setGraph(new PipeViewGraph(getGraphModel()));
    //getView().connectToMouseLeftClickRouter(getGraph(), getView().SELECT_GRAPH_NODE_PARENTS_KEY);
    GridBagConstraints constraints = new GridBagConstraints();
    constraints.gridx = 0;
    constraints.gridy = 0;
    constraints.weightx = 1;
    constraints.weighty = 1;
    constraints.gridwidth = 3;
    constraints.fill = GridBagConstraints.BOTH;
    setScrollPane(new JScrollPane(getGraph()));
    add(getScrollPane(), constraints);
    getGraph().setEditable(false);
    getGraph().setMoveable(true);
    getGraph().setCloneable(false);
    getGraph().setBendable(false);
    getGraph().setConnectable(false);
    getGraph().setDisconnectable(false);
    getView().reshow();
  }
  
  private void storeKeysAndPositions(java.util.List list, int row){
    ModelElement element = null;
    for(int column=0; column<list.size(); column++){
      element = (ModelElement)list.get(column);
      if(!getKeyToPositionMap().containsKey(element.getKey())){
        if(getView().getLogger().isLoggingHigh()){
          getView().getLogger().logHigh("key: "+element.getKey()+" storing position: "+row+" "+column);
        }
        getKeyToPositionMap().put(
          element.getKey(), 
          new Dyad(column, row)
        );
      }else{
        if(getView().getLogger().isLoggingHigh()){
          getView().getLogger().logHigh("key: "+element.getKey()+" already appears - not storing");
        }
      }
    }
  }
  
  private ArrayList listChildren(Collection elements){
    HashSet returnSet = new HashSet();
    ArrayList returnList = new ArrayList();
    Iterator iterator = elements.iterator();
    ModelElement element = null;
    while(iterator.hasNext()){
      element = (ModelElement)iterator.next();
      returnSet.addAll(element.getChildElements().values());
    }
    returnList.addAll(returnSet);
    return returnList;
  }

  /**
   * Draw a row of rectangles, one for each element of the list, at the
   * indicated row
  **/
  private ArrayList recordEdges(int iterate, java.util.List list){
    Iterator elements = list.iterator();
    ModelElement element;
    Iterator parents;
    ModelElement parent;
    ArrayList returnList = new ArrayList();
    Link link = null;
    while(elements.hasNext()){
      element = (ModelElement)elements.next();
      if(iterate > 0){
        parents = element.getParentElements().iterator();
        while(parents.hasNext()){
          parent = (ModelElement)parents.next();
          link = new Link();
          link.from = parent.getKey();
          link.to = element.getKey();
          returnList.add(link);
          if(getView().getLogger().isLoggingHigh()){
            getView().getLogger().logHigh("Created graph edge between "+parent.getKey()+" and "+element.getKey());
          }//end if
        }//end while
      }//end if
    }//end for
    
    return returnList;
  }//end recordEdges
  
  private void drawLinks(int iterate, java.util.List list){
    Iterator elements = list.iterator();
    ModelElement element;
    Iterator parents;
    ModelElement parent;
    while(elements.hasNext()){
      element = (ModelElement)elements.next();
      if(iterate > 0){
        parents = element.getParentElements().iterator();
        while(parents.hasNext()){
          parent = (ModelElement)parents.next();
          drawLineBetween(parent, element);
          if(getView().getLogger().isLoggingHigh()){
            getView().getLogger().logHigh("Created graph edge between "+parent.getKey()+" and "+element.getKey());
          }//end if
        }//end while
      }//end if
    }//end for
  }
  
  private void drawLineBetween(ModelElement parent, ModelElement child){
    DefaultEdge edge = new DefaultEdge();
    Port source = (Port)getElementPortMap().get(parent.getKey());
    Port target = (Port)getElementPortMap().get(child.getKey());
    ConnectionSet connectionSet = new ConnectionSet(edge, source, target);
    Map edgeAttributes = GraphConstants.createMap();
    Object[] insert = new Object[]{edge};
    HashMap tempMap = new HashMap();
    int arrow = GraphConstants.ARROW_CLASSIC;
    GraphConstants.setLineEnd(edgeAttributes, arrow);
    GraphConstants.setEndFill(edgeAttributes, true);
    tempMap.put(edge, edgeAttributes);
    
    getGraphModel().insert(insert, tempMap, connectionSet, null, null);
    
  }
  
  private void drawRectangle(
    ModelElement element,
    int x, 
    int y, 
    int width, 
    int height, 
    String elementKey
  ){
    DefaultGraphCell cell = new DefaultGraphCell(elementKey);
    Map cellMap = GraphConstants.createMap();
    GraphConstants.setBounds(cellMap, new Rectangle(x, y, width,height));
    GraphConstants.setBorderColor(cellMap, Color.black);
    GraphConstants.setOpaque(cellMap, true);
    
    //copy over all the properties of this cell node in the model into 
    //the properties of the graph model's node.
    cellMap.putAll(element.getProperties());
    
    DefaultPort port;
    HashMap cellToAttributesMap = new HashMap();
      
    cellToAttributesMap.put(cell, cellMap);
    port = new DefaultPort();
    cell.add(port);
    Object[] insert = new Object[]{cell};
    getGraphModel().insert(insert,  cellToAttributesMap, null, null, null);
    getElementPortMap().put(elementKey, port);
    getGraphAttributeMap().put(elementKey, cellMap);
    getCellMap().put(elementKey, cell);
  }
  
  public void update(PathStepsModel model){
    ModelElement rootElement = model.getRootElement();
    ModelElement panelModel = rootElement.getChildElement(model.PATH_STEPS_PANEL);
    
    if(panelModel == null){
      return;
    }
    
    Collection allNodes = (Collection)panelModel.getProperty(model.PATH_STEPS_PANEL_ALL_NODES);
    Iterator nodes = allNodes.iterator();
    ModelElement node;
    Map cellAttributes;
    Rectangle cellBounds;
    int cellX = 0;
    int cellY = 0;
    int cellWidth = 0;
    int cellHeight = 0;
    String boundString = null;
    DefaultGraphCell cell = null;
    
    Properties graphLayout = 
      (Properties)panelModel
        .getProperty(
          model.PATH_STEPS_PANEL_GRAPH_LAYOUT_CONFIGURATION
        );
    
    if(graphLayout == null){
      throw new FatalAException(
        "Expected a model (with a Properties value) of the element property: "+
        model.PATH_STEPS_PANEL_GRAPH_LAYOUT_CONFIGURATION+
        " for the model element :"+model.PATH_STEPS_PANEL_ALL_NODES
      );  
    }

    while(nodes.hasNext()){
      node = (ModelElement)nodes.next();
      cellAttributes = (Map)((DefaultGraphCell)getCellMap().get(node.getKey())).getAttributes();
      //cellAttributes = (Map)getGraphAttributeMap().get(node.getKey());
      if(cellAttributes == null){
        throw new FatalAException("Cannot find the attribute map for cell of key "+node.getKey());
      }
      cellBounds = GraphConstants.getBounds(cellAttributes); 
      cellX = new Double(cellBounds.getX()).intValue();
      cellY = new Double(cellBounds.getY()).intValue();
      cellWidth = new Double(cellBounds.getWidth()).intValue();
      cellHeight = new Double(cellBounds.getHeight()).intValue();
      boundString = 
        (new StringBuffer())
          .append(String.valueOf(cellX))
          .append(" ")
          .append(String.valueOf(cellY))
          .append(" ")
          .append(String.valueOf(cellWidth))
          .append(" ")
          .append(String.valueOf(cellHeight))
          .toString();
       graphLayout.setProperty(node.getKey(), boundString);
    }
  }
  
  public void setModel(PathStepsModel model){
    this.model = model;
  }
  
  private PathStepsModel getModel(){
    return model;
  }
  
  private HashMap getKeyToPositionMap(){
    return keyToPositionMap;
  }
  
  private PathStepsSwingView getView(){
    return view;
  }
  
  private void setGraphModel(GraphModel model){
    this.graphModel = model;
  }
  
  private GraphModel getGraphModel(){
    return graphModel;
  }
  
  private void setGraph(JGraph graph){
    this.graph = graph;
  }
  
  private JGraph getGraph(){
    return graph;
  }
  
  private HashMap getGraphAttributeMap(){
    return graphAttributeMap;
  }
  
  private HashMap getCellMap(){
    return cellMap;
  }
  
  private HashMap getElementPortMap(){
    return elementPortMap;
  }
  
  private FontMetrics getFontMetrics(){
    if(fontMetrics == null){
      fontMetrics = getGraphics().getFontMetrics();
    }
      
    return fontMetrics;
  }
  
  public void clear(){
    if(getView().getLogger().isLoggingHigh()){
      getView().getLogger().logHigh("Starting clear");
    }
    remove(scrollPane);
    if(getView().getLogger().isLoggingHigh()){
      getView().getLogger().logHigh("Finishing clear");
    }
  }
  
  private JScrollPane getScrollPane(){
    return scrollPane;
  }
  
  private void setScrollPane(JScrollPane pane){
    scrollPane = pane;
  }
  
  private HashSet getDrawnSet(){
    return drawnSet;
  }
  
  private boolean getBooleanFromModel(ModelElement element, String key){
    String booleanString = (String)element.getProperty(key);

    boolean returnValue;
    
    if(booleanString == null || booleanString.trim().length() <=0){
      throw new FatalAException("You must provide a value in the application config for: "+key);
    }

    if(
      !Boolean.FALSE.toString().equals(booleanString) &&
      !Boolean.TRUE.toString().equals(booleanString)
    ){
      throw new FatalAException("The value "+key+" could not be parsed into a boolean");
    }
    
    returnValue = Boolean.valueOf(booleanString).booleanValue();
    
    return returnValue;
  }
  
  private int getNumberFromModel(ModelElement element, String key){
    String numberString = (String)element.getProperty(key);

    int number = 0;
    if(numberString == null || numberString.trim().length() <=0){
      throw new FatalAException("You must provide a value in the application config for: "+key);
    }

    try{
      number = Integer.valueOf(numberString).intValue();
    }catch(NumberFormatException exception){
      throw new FatalAException(
        "The layout config item: "+key+" is "+numberString+
        " which is not a valid integer", 
        exception
      );
    }
    
    return number;
  }

  public Object getSelectedGraphCellObject(){
    if(getGraph() == null){
      throw new FatalAException("Attempt made to access graph in PathStepsPanel before it was initialsed");
    }
    return getGraph().getSelectionCell();
  }
  
  public void applyGraphLayout(){

    ArrayList list = new ArrayList();
    Set keySet = getKeyToPositionMap().keySet();
    Iterator keys = keySet.iterator();
    String key = null;
    Node node = null;
    Dyad position = null;
    ArrayList rootKeys = new ArrayList();
    ArrayList allEdges = new ArrayList();
    LayoutCalculator calculator = null;
    GraphUpdater updater;
    ModelElement rootElement = null;
    ModelElement layoutDialogElement = null;
    int row = 0;
    Object[] roots;
    Comparator inverseChildNumber =  new PathStepsPanel.InverseChildNumberComparator();
    
    int xSeparation = 0;
    int ySeparation = 0;
    int naturalEdgeLength = 0;
    int repulsionMultiplier = 0;
    int iterateNumber = 0;
    int movementLimit = 0;
    int gravity = 0;
    boolean fixRoots = false;

    Map cellAttributes = null;
    Rectangle cellBounds = null;
      
    if(getView().getLogger().isLoggingHigh()){
      getView().getLogger().logHigh("Started applying graph layout");
    }
    
    if(model == null){
      throw new FatalAException("Attempt to apply graph layout to a null model");
    }
    rootElement = model.getRootElement().getChildElement(PathStepsModel.PATH_STEPS_PANEL);
    if(rootElement == null){
      throw new FatalAException("Attempt to apply graph layout when PATH_STEPS_PANEL root is null");
    }
    
    layoutDialogElement = model.getRootElement().getChildElement(PathStepsModel.LAYOUT_DIALOG);
    xSeparation = getNumberFromModel(layoutDialogElement, model.LAYOUT_DIALOG_HORIZONTAL_SPACING);
    ySeparation = getNumberFromModel(layoutDialogElement, model.LAYOUT_DIALOG_VERTICAL_SPACING);
    naturalEdgeLength = getNumberFromModel(layoutDialogElement, model.LAYOUT_DIALOG_SPRING_NATURAL_LENGTH);
    repulsionMultiplier = getNumberFromModel(layoutDialogElement, model.LAYOUT_DIALOG_REPULSION_MULTIPLIER);
    iterateNumber = getNumberFromModel(layoutDialogElement, model.LAYOUT_DIALOG_ITERATES);
    movementLimit = getNumberFromModel(layoutDialogElement, model.LAYOUT_DIALOG_MOVEMENT_LIMIT);
    gravity = getNumberFromModel(layoutDialogElement, model.LAYOUT_DIALOG_GRAVITY);
    fixRoots = getBooleanFromModel(layoutDialogElement, model.LAYOUT_DIALOG_FIX_ROOTS);

    list.clear();
    list.add(rootElement);
    list = listChildren(list);

    //list contains the first row - the 'roots'. store their keys.
    for(int i=0; i<list.size(); i++){
      rootKeys.add(((ModelElement)list.get(i)).getKey());
    }

    int counter=0;
    nodes = new Node[keySet.size()];
    while(keys.hasNext()){
      key = (String)keys.next();
      position = (Dyad)getKeyToPositionMap().get(key);
      node = new Node();


      cellAttributes = ((DefaultGraphCell)getCellMap().get(key)).getAttributes();
      if(cellAttributes == null){
        throw new FatalAException("Cannot find the attribute map for cell of key "+key);
      }

      cellBounds = GraphConstants.getBounds(cellAttributes);
      System.out.println(
        " <<<< retrieved key "+key+" with cell: "+getCellMap().get(key).hashCode()+" with bounds: "+
        GraphConstants.getBounds(((DefaultGraphCell)getCellMap().get(key)).getAttributes())
      );
      
      node.x = cellBounds.x;
      node.y = cellBounds.y;
      node.lbl = key;
      
      if(getView().getLogger().isLoggingHigh()){
        getView().getLogger().logHigh("node: "+key+":"+node.x+"-"+node.y);
      }

      //make fixed if root.
      if(fixRoots && rootKeys.contains(key)){
        if(getView().getLogger().isLoggingHigh()){
          getView().getLogger().logHigh("FIXING node with key: "+key);
        }
        node.fixed = true;
      }
      nodes[counter]=node;
      counter++;
    }

    //
    //Gather up all edges
    row = 0;
    list.clear();
    list.add(rootElement);
    list = listChildren(list);
    Collections.sort(list, inverseChildNumber);
    while(!list.isEmpty()){
      allEdges.addAll(recordEdges(row, list));
      Collections.sort(list, inverseChildNumber);
      list = listChildren(list);
      row++;
    }

    //
    //Go through each edge and add to the Edge[] array
    pathsteps.gui.Edge edge = null;
    Link link = null;
    edges = new pathsteps.gui.Edge[allEdges.size()];
    for(int i=0; i<allEdges.size(); i++){
      link = (Link)allEdges.get(i);
      edge = new pathsteps.gui.Edge();
      edge.from = findNode(link.from);
      edge.to = findNode(link.to);
      edge.len = naturalEdgeLength;
      edges[i] = edge;
    }
    
    calculator = new LayoutCalculator(nodes, edges, rootKeys.size() * xSeparation, 1000);

    updater = 
      new GraphUpdater(
        calculator, 
        getGraphModel(),
        getCellMap(),
        getGraphAttributeMap(),
        iterateNumber, 
        1, 
        movementLimit, 
        gravity, 
        repulsionMultiplier
      );

    javax.swing.Timer timer = new javax.swing.Timer(
      10, 
      updater
    );

    updater.timer = timer;

    if(getView().getLogger().isLoggingMedium()){
      getView().getLogger().logMedium("Starting timer");
    }
    timer.start();
  }
  
  private Dyad parsePositionString(String position){
    StringTokenizer tokenizer = new StringTokenizer(position);
    String number;
    int x = 0;
    int y = 0;
    int count = 0;
    while(tokenizer.hasMoreTokens() && count < 2){
      number = tokenizer.nextToken();
      if(count == 0){
        x = Integer.valueOf(number).intValue();
      }else{
        y = Integer.valueOf(number).intValue();
      }
      count++;
    }
    
    if(count < 1){
      throw new FatalAException("Found a graph layout configuration with at less than two coordinates per node!");
    }
    return new Dyad(x,y);
  }
}