package pathsteps.gui;
import pathsteps.common.*;
import pathsteps.model.*;
import java.util.*;

public class SetUtil {
  ALogger _logger = null;
  Comparator degreeComparator = new DegreeComparator();
  
  public class DegreeComparator implements Comparator{
    public int compare(Object node1, Object node2) {
      ModelElement element1 = (ModelElement)node1;
      ModelElement element2 = (ModelElement)node2;
      Object object1 = (Integer)element1.getProperty(ModelElement.DEGREE);
      Object object2 = (Integer)element2.getProperty(ModelElement.DEGREE);
      
      if(!(object1 instanceof Integer)|| !(object2 instanceof Integer)){
        throw new FatalAException("DEGREE property of ModelElement are "+object1+" / "+object2+" - not Integer in degree comparator");
      }
      
      Integer degree1 = (Integer)object1;
      Integer degree2 = (Integer)object2;
      
      if(degree1 == null || degree2 == null){
        throw new FatalAException("Attempt to access non-existent DEGREE property in degree comparator");
      }
      
      return degree1.compareTo(degree2);
    }
    
    public boolean equals(Object node1, Object node2){
      ModelElement element1 = (ModelElement)node1;
      ModelElement element2 = (ModelElement)node2;
      Object object1 = (Integer)element1.getProperty(ModelElement.DEGREE);
      Object object2 = (Integer)element2.getProperty(ModelElement.DEGREE);
      
      if(!(object1 instanceof Integer)|| !(object2 instanceof Integer)){
        throw new FatalAException("DEGREE property of ModelElement are "+object1+" / "+object2+" - not Integer in degree comparator");
      }
      
      Integer degree1 = (Integer)object1;
      Integer degree2 = (Integer)object2;
      
      if(degree1 == null || degree2 == null){
        throw new FatalAException("Attempt to access non-existent DEGREE property in degree comparator");
      }
      
      return degree1.equals(degree2);
    }    
  }
  
  public SetUtil(ALogger logger) {
    _logger = logger;
  }

  /**
   * Assumes you've been passed a connected graph, and assigns x-y positions to elements
   * of graph based on a depth-first search
  **/
  public void positionElementsOfSet(Set returnSet, ModelElement seed, int xDegree, int yDegree){
    ModelElement relatedElement;
    Collection childList;

    if(!returnSet.contains(seed)){
      seed.addProperty(ModelElement.X, new Integer(xDegree));
      seed.addProperty(ModelElement.Y, new Integer(yDegree));
      returnSet.add(seed);
      if(getLogger().isLoggingMedium()){
        getLogger().logMedium("Positioning element "+seed.getKey()+" at "+xDegree+", "+yDegree);
      }

      childList = (new ArrayList());
      childList.addAll(seed.getChildElements().values());
      Iterator kids = childList.iterator();
      while(kids.hasNext()){
        relatedElement = (ModelElement)kids.next();
        positionElementsOfSet(returnSet, relatedElement, xDegree, yDegree+1);
      }

      childList.clear();
      childList.addAll(seed.getChildElements().values());
      kids = childList.iterator();
      int childCounter = 0;
      while(kids.hasNext()){
        relatedElement = (ModelElement)kids.next();
        if(!relatedElement.getKey().equals(PathStepsModel.PATH_STEPS_PANEL)){
          positionElementsOfSet(returnSet, relatedElement, xDegree + 1, yDegree);
          childCounter++;
        }
      }
    }
  }
  
  public Dyad getMaxExtentsOfSet(Set set){
    int x = 0;
    int y = 0;
    int maxX = 0;
    int maxY = 0;
    ModelElement element;
    Iterator elements = set.iterator();
    while(elements.hasNext()){
      element = (ModelElement)elements.next();
      x = ((Integer)(element.getProperty(ModelElement.X))).intValue();
      y = ((Integer)(element.getProperty(ModelElement.X))).intValue();
      if(x > maxX){
        maxX = x;
      }
      
      if(y > maxY){
        maxY = y;
      }
    }
    return new Dyad(x,y);
  }
  
  public List positionElementOfList(List modelElementList){
    ModelElement element;
    Iterator elements = modelElementList.iterator();
    List returnList = new ArrayList();
    Set currentSet = null;
    Set allFoundElements = new HashSet();
    Dyad position;
    int x = 0, y = 0;
    
    while(elements.hasNext()){
      element = (ModelElement)elements.next();
      if(getLogger().isLoggingMedium()){
        getLogger().logMedium("Positioning elements of list: considering element "+element.getKey());
      }
      
      if(!allFoundElements.contains(element)){
        if(currentSet != null){
          position = getMaxExtentsOfSet(currentSet);
          x = position.x();
          y = position.y();
          if(getLogger().isLoggingMedium()){
            getLogger().logMedium("x-y position of last positioned set is: "+x+", "+y);
          }
        }
        
        currentSet = new HashSet();
        if(getLogger().isLoggingMedium()){
          getLogger().logMedium("element not found: creating set "+currentSet.hashCode());
        }
        
        positionElementsOfSet(currentSet, element, x + 1, 0);
        allFoundElements.addAll(currentSet);
        returnList.add(currentSet);
      }else{
        if(getLogger().isLoggingMedium()){
          getLogger().logMedium("element has already been found");
        }
      }
    }
    
    return returnList;
  }
  
  /**
   * Returns the connected components of an input graph, as a List of Sets
   **/
  public void findAllConnectedElements(Set returnSet, ModelElement seed, boolean markDegree, int degree){
    if(!returnSet.contains(seed)){
      ModelElement relatedElement;
      if(markDegree){
        seed.addProperty(ModelElement.DEGREE, new Integer(degree));
      }
      returnSet.add(seed);
      if(getLogger().isLoggingMedium()){
        getLogger().logMedium("Adding element "+seed.getKey()+" to set "+returnSet.hashCode());
      }
      Collection childList = (new ArrayList());
      childList.addAll(seed.getChildElements().values());
      childList.addAll(seed.getParentElements());
      
      Iterator kids = childList.iterator();
      while(kids.hasNext()){
        relatedElement = (ModelElement)kids.next();
        if(!relatedElement.getKey().equals(PathStepsModel.PATH_STEPS_PANEL)){
          findAllConnectedElements(returnSet, relatedElement, markDegree, degree+1);
        }
      }
    }
  }

  /**
   * Input - a list of ModelElements (presumably the root of some tree).
   * Output - a List of Sets, each containing the connected ModelElements
   **/
  public List findConnectedComponents(List modelElementList){
    ModelElement element;
    Iterator elements = modelElementList.iterator();
    List returnList = new ArrayList();
    Set currentSet;
    Set allFoundElements = new HashSet();
    
    while(elements.hasNext()){
      element = (ModelElement)elements.next();
      if(getLogger().isLoggingMedium()){
        getLogger().logMedium("Considering element "+element.getKey());
      }
      if(!allFoundElements.contains(element)){
        currentSet = new HashSet();
        if(getLogger().isLoggingMedium()){
          getLogger().logMedium("element not found: creating set "+currentSet.hashCode());
        }
        findAllConnectedElements(currentSet, element,true, 0);
        allFoundElements.addAll(currentSet);
        returnList.add(currentSet);
      }else{
        if(getLogger().isLoggingMedium()){
          getLogger().logMedium("element has already been found");
        }
      }
    }
    
    return returnList;
  }
  
  public List sortRootsByConnectivity(List roots){
    List componentList = findConnectedComponents(roots);
    List returnList = new ArrayList();
    List setList = new ArrayList();
    Set set;
    Iterator elements;
    Object element;
    for(int setCounter=0; setCounter<componentList.size(); setCounter++){
      set = (Set)componentList.get(setCounter);
      setList.clear();
      elements = roots.iterator();
      while(elements.hasNext()){
        element = elements.next();
        if(set.contains(element)){
          setList.add(element);
        }
      }
      Collections.sort(setList, getDegreeComparator());
      returnList.addAll(setList);
    }
    return returnList;
  }
  
  private ALogger getLogger(){
    return _logger;
  }
  
  public Comparator getDegreeComparator(){
    return degreeComparator;
  }
  
  public void clearPropertyFromChildrenOfNode(ModelElement root, String propertyKey){
    Set allElements = new HashSet();
    //
    //Collect all elements into the allElements set
    findAllConnectedElements(allElements, root, false, 0);
    //
    //Now walk the elements and clear away the SELECTED property
    Iterator elements = allElements.iterator();
    while(elements.hasNext()){
      ((ModelElement)elements.next()).removeProperty(propertyKey);
    }
  }
  
  public void addPropertyToParentsOfNode(ModelElement startElement, ModelElement rootOfTree, String key, String value){
    Set parentSet = new HashSet();
    ModelElement parentElement;
    parentSet.add(startElement);
    startElement.addProperty(key, value);
    Iterator parents = startElement.getParentElements().iterator();
    while(parents.hasNext()){
      parentElement = (ModelElement)parents.next();
      if(!(parentElement == rootOfTree)){
        addPropertyToParentsOfNode((ModelElement)parentElement, rootOfTree, key, value);
      }
    }
  }
  
  public Set findAllChildren(ModelElement root){
    Set returnSet = new HashSet();
    return returnSet;
  }
  
  public ModelElement findRelatedElementWithName(ModelElement treeRoot, String name){
    Set allElements = new HashSet();
    ModelElement element;
    //
    //Collect all elements into the allElements set
    findAllConnectedElements(allElements, treeRoot, false, 0);

    Iterator elements = allElements.iterator();
    while(elements.hasNext()){
      element = (ModelElement)elements.next();
      if(element.getKey().equals(name)){
        return element;
      }
    }
    
    return null;
  }
  
  public static void main(String[] args){
    ModelElement rootElement = new ModelElement(PathStepsModel.PATH_STEPS_PANEL);
    ModelElement e1 = rootElement.createChildElement("e1");
    ModelElement e2 = rootElement.createChildElement("e2");
    ModelElement e3 = rootElement.createChildElement("e3");
    ModelElement e4 = rootElement.createChildElement("e4");
    
    ModelElement e11 = e1.createChildElement("e11");
    ModelElement e111 = e11.createChildElement("e111");
    
    ModelElement e31 = e3.createChildElement("e31");
    ModelElement e32 = e3.createChildElement("e32");
    
    ModelElement e41 = e4.createChildElement("e41");
    ModelElement e42 = e4.createChildElement("e42");
    
    e111.addParentElement(e2);
    
    e1.addChildElement(e42);
    
    List input = new ArrayList();
    HashSet inputSet = new HashSet();
    
    input.add(e1);
    input.add(e3);
    input.add(e2);
    input.add(e4);
    
    inputSet.add(e1);
    inputSet.add(e3);
    inputSet.add(e2);
    inputSet.add(e4);
    
    SetUtil util = new SetUtil(new SimpleLogger());
    /*
    List returnList = util.findConnectedComponents(input);
    for(int i=0; i<returnList.size();i++){
      System.out.println("Set size: "+((Set)returnList.get(i)).size());
      Set set = (Set)returnList.get(i);
      Iterator iterator = set.iterator();
      while(iterator.hasNext()){
        System.out.println(iterator.next());
      }
    }
    
    List sortedRoots = util.sortRootsByConnectivity(input);
    for(int i=0; i<sortedRoots.size();i++){
      System.out.println("Sorted root: "+((ModelElement)sortedRoots.get(i)).getKey());
    }
     */
    
    List list = util.positionElementOfList(input);
    Set set = null;
    /*
    for(int i=0; i<list.size();i++){
      System.out.println("Sorted set: "+((Set)list.get(i)).hashCode());
      Iterator elements = set.iterator();
      while(elements.hasNext()){
        System.out.println("Element: "+((ModelElement)elements.next()).getKey());
      }
    }
     */
  }
  
}
