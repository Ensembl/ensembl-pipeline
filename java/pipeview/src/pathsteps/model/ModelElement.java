package pathsteps.model;
import java.util.*;

import pathsteps.common.*;

/**
 * A generic composite class to build a model from. There are two types of children 
 * contained: key-value string pairings (properties) and child ModelElements (also 
 * stored as key-value pairs).
**/
public class ModelElement{
  private String key;
  private HashMap parentElements = new HashMap();

  private HashMap childElements = new HashMap();
  private HashMap properties = new HashMap();
  private boolean dirty;
  private boolean _new = true;
  public static String DEGREE = "DEGREE";
  public static String SELECTED = "SELECTED";
  public static String X = "X";
  public static String Y = "Y";
  
  public static String ANALYSIS_STATUS_COUNT = "ANALYSIS_STATUS_COUNT";
  public static String FINISHED_COUNT = "FINISHED_COUNT";
  public static String CREATED = "CREATED";
  public static String FAILED = "FAILED";
  public static String READING = "READING";
  public static String RUNNING = "RUNNING";
  public static String SUBMITTED = "SUBMITTED";
  public static String SHOW_DETAIL = "SHOW_DETAIL";
  
  public ModelElement(String key){
    this.key = key;
  }
  
  public ModelElement(ModelElement parent, String key){
    addParentElement(parent);
    this.key = key;
  }
  
  public String getKey(){
    return key;
  }
  
  public HashMap getChildElements(){
    return childElements;
  }
  
  public HashMap getProperties(){
    return properties;
  }
  
  public Collection getPropertyValues(){
    return getProperties().values();
  }
  
  public Set getPropertyKeys(){
    return getProperties().keySet();
  }
  
  public void addChildElement(String key, ModelElement element){
    getChildElements().put(key,  element);
    element.addParentElement(this);
  }
  
  public void addChildElement(ModelElement element){
    getChildElements().put(element.getKey(),  element);
    element.addParentElement(this);
  }

  public void removeChildElement(String key){
    ModelElement element = getChildElement(key);
    if(element == null){
      throw new FatalAException("Attempt to remove non-existent child element "+key+" from parent "+getKey());
    }
    getChildElements().remove(key);
  }
  
  public ModelElement createChildElement(String key){
    ModelElement element =  new ModelElement(this, key);
    getChildElements().put(key, element);
    return element;
  }
  
  public void addProperty(String key, Object value){
    getProperties().put(key, value);
  }
  
  public void removeProperty(String key){
    getProperties().remove(key);
  }
  
  public ModelElement getChildElement(String key){
    return (ModelElement)getChildElements().get(key);
  }
  
  public Object getProperty(String key){
    return getProperties().get(key);
  }
  
  public void setDirty(boolean value){
    dirty = value;
  }
  
  public void setNew(boolean value){
    _new = value;
  }
  
  public void dirty(){
    setDirty(true);
  }
  
  public void clean(){
    setDirty(false);
  }
  
  public boolean isDirty(){
    return dirty;
  }
  
  public boolean isNew(){
    return _new;
  }
  
  public void renew(){
    setNew(true);
  }
  
  public void age(){
    setNew(false);
  }
  
  public String toString(){
    return toString(0);
  }

  public void addParentElement(ModelElement parent){
    privateGetParentElements().put(parent.getKey(), parent);
  }
  
  public Collection getParentElements(){
    return privateGetParentElements().values();
  }
  
  public ModelElement getParentElement(String key){
    return (ModelElement)privateGetParentElements().get(key);
  }
  
  private HashMap privateGetParentElements(){
    return parentElements;
  }
  
  public String getParentString(){
    Iterator parents = getParentElements().iterator();
    StringBuffer buffer = new StringBuffer().append("(");
    while(parents.hasNext()){
      buffer.append(((ModelElement)parents.next()).getKey()).append(",");
    }
    buffer.append(")");
    return buffer.toString();
  }  

  public String toString(int depth){
    StringBuffer buffer = new StringBuffer();
    String key = null;
    Object value = null;
    String indent="";
    Iterator keys = null;
    String retrn = "\n";
    for(int i=0;i<depth;i++){
      indent+="    ";
    }
    
    buffer
      .append(indent)
      .append(getKey())
      .append("(");
    
    if(getParentElements().size()>0){
      buffer.append(getParentString());
    }
    
    buffer.append(retrn);
    
    keys = getProperties().keySet().iterator();
    while(keys.hasNext()){
      key = (String)keys.next();
      value = getProperty(key);
      buffer.append(indent).append(key).append(":").append(value).append(retrn);
    }
    
    keys = getChildElements().keySet().iterator();
    while(keys.hasNext()){
      key = (String)keys.next();
      value = getChildElement(key);
      buffer.append(((ModelElement)value).toString(depth+1));
    }
    
    return buffer.toString();
  }
}
