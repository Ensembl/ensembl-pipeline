package pathsteps.model;
import java.util.*;

import pathsteps.common.*;

/**
 * An extension of the ModelElement class which allows an element to have multiple parents.
**/
public class ExtendedModelElement
extends ModelElement
{
  private HashMap parentElements = new HashMap();
  
  public ExtendedModelElement(String key){
    super(key);
  }
  
  public ExtendedModelElement(ModelElement parent, String key){
    super(parent, key);
  }
}
