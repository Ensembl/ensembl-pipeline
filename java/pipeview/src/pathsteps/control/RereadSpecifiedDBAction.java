package pathsteps.control;

import pathsteps.common.*;
import pathsteps.model.*;

import java.util.*;
import java.sql.*;

/**
 * When the user types into the db-change field, the databases in the dropdown
 * should be blanked out: this is what this action does.
**/
public class RereadSpecifiedDBAction extends AAction{

  public RereadSpecifiedDBAction(AEventHandler eventHandler){
    super(eventHandler);
  }
  
  protected void doAction(AView view, AModel aModel){
    PathStepsModel model = (PathStepsModel)aModel;
    String host = null;
    String port = null;
    String user = null;
    String password = null;
    Collection databases = null;
    String selectedDatabase = null;
    ModelElement dialogModel = model.getRootElement().getChildElement(model.READ_DB_DIALOG);
    Collection list;
    
    if(getLogger().isLoggingLow()){
      getLogger().logLow("Starting ReadSpecifiedDBAction");
    }
    
    if(dialogModel == null){
      throw new FatalAException("ReadDBDialog expecting a model called: "+model.READ_DB_DIALOG+" but none was found");
    }
    
    host = (String)dialogModel.getProperty(model.READ_DB_DIALOG_HOST);
    port = (String)dialogModel.getProperty(model.READ_DB_DIALOG_PORT);
    user = (String)dialogModel.getProperty(model.READ_DB_DIALOG_USER);

    password = (String)dialogModel.getProperty(model.READ_DB_DIALOG_PASSWORD);
    selectedDatabase = (String)dialogModel.getProperty(model.READ_DB_DIALOG_PIPELINE_DB);
    
    if(user == null || user.trim().length()<=0){
      throw new NonFatalAException("User must be supplied");
    }
    
    if(host == null || host.trim().length()<=0){
      throw new NonFatalAException("Host must be supplied");
    }
    
    if(port == null || port.trim().length()<=0){
      throw new NonFatalAException("Port must be supplied");
    }
    
    if(selectedDatabase == null || selectedDatabase.trim().length()<=0){
      throw new NonFatalAException("Ensembl Database must be supplied");
    }
    
    //check that the chosen db has all the tables we want, throw an exception if not.
    //Run the logic to store the chains of dependent path/steps in the model
    model.getRootElement().addChildElement(
      PathStepsModel.PATH_STEPS_PANEL,
      createPathStepsModel(selectedDatabase)
    );
    
    //System.out.println(model.toString());
    Properties history = view.getApplication().getHistory();
    history.setProperty(model.READ_DB_DIALOG_HOST, host);
    history.setProperty(model.READ_DB_DIALOG_PORT, port);
    history.setProperty(model.READ_DB_DIALOG_USER, user);
    
    getView().getApplication().writeHistory(history);
    
    if(getLogger().isLoggingLow()){
      getLogger().logLow("Finished ReadSpecifiedDBAction");
    }
    
    if(getLogger().isLoggingLow()){
      getLogger().logLow("Wrote application history");
    }
  }

  private ModelElement createFakePathStepsModel(String selectedDatabase){
    ModelElement rootElement = new ModelElement(PathStepsModel.PATH_STEPS_PANEL);
    /*
    ModelElement submitContig = rootElement.createChildElement("SubmitContig");
    ModelElement repeatMask = submitContig.createChildElement("RepeatMask");
    ModelElement CpG = submitContig.createChildElement("CpG");
    ExtendedModelElement accumulator = new ExtendedModelElement("Accumulator");
    repeatMask.addChildElement(accumulator);
    CpG.addChildElement(accumulator);
     */
    ModelElement e1 = rootElement.createChildElement("e1");
    ModelElement e2 = rootElement.createChildElement("e2");
    ModelElement e3 = rootElement.createChildElement("e3");
    ModelElement e4 = rootElement.createChildElement("e4");
    
    ModelElement e11 = e1.createChildElement("e11");
    e11.addParentElement(e2);
    
    ModelElement e111 = e11.createChildElement("e111");
    
    ModelElement e1111 = e111.createChildElement("e1111");

    ModelElement e31 = e3.createChildElement("e31");
    ModelElement e32 = e3.createChildElement("e32");
    
    ModelElement e41 = e4.createChildElement("e41");
    ModelElement e42 = e4.createChildElement("e42");    
    e111.addParentElement(e41);

    if(getLogger().isLoggingMedium()){
      getLogger().logMedium("Finished creating model");
    }
    
    List elements = new ArrayList();
    elements.add(e1);
    elements.add(e2);
    elements.add(e3);
    elements.add(e11);
    elements.add(e111);
    elements.add(e3);
    elements.add(e31);
    elements.add(e32);
    elements.add(e4);
    elements.add(e41);
    elements.add(e42);
    
    rootElement.addProperty(PathStepsModel.PATH_STEPS_PANEL_ALL_NODES, elements);
    
    return rootElement;
  }
  
  private ModelElement createPathStepsModel(String selectedDatabase){
    
    ModelElement rootElement = new ModelElement(PathStepsModel.PATH_STEPS_PANEL);
    Statement statement;
    ResultSet resultSet;
    String select;
    HashMap elementsByAnalysisId = new HashMap();
    HashMap elementsByLogicName = new HashMap();
    HashMap elementsByRuleId = new HashMap();
    ExtendedModelElement currentElement;
    int analysis_id;
    String logic_name;
    int rule_id;
    String rule_condition;
    ModelElement parentElement;
    
    if(getLogger().isLoggingMedium()){
      getLogger().logMedium("Started creating model");
    }

    
    try{
      java.sql.Connection conn = createNewConnectionFromSelectedDatabase((PathStepsModel)getModel());
      statement = conn.createStatement();
      
      //
      //Read in all analyses - create a unique element for each analysis (unique
      //by logic name)
      select = "select analysis_id, logic_name from analysis";
      resultSet = statement.executeQuery(select);
      while(resultSet.next()){
        analysis_id = resultSet.getInt("analysis_id");
        logic_name = resultSet.getString("logic_name");
        currentElement = new ExtendedModelElement(logic_name);
        currentElement.addProperty("analysis_id", String.valueOf(analysis_id));
        currentElement.addProperty("logic_name", String.valueOf(logic_name));
        elementsByLogicName.put(logic_name, currentElement);
        elementsByAnalysisId.put(String.valueOf(analysis_id), currentElement);
      }
      
      //
      //Read in the rules which point at these analyses
      select = "select rule_id, goal from rule_goal";
      resultSet = statement.executeQuery(select);
      while(resultSet.next()){
        rule_id = resultSet.getInt("rule_id");
        analysis_id = resultSet.getInt("goal");
        currentElement = (ExtendedModelElement)elementsByAnalysisId.get(String.valueOf(analysis_id));
        if(currentElement == null){
          throw new FatalAException("Rule Goal: "+rule_id+" refers to analysis "+analysis_id+" which doesn't exist");
        }
        elementsByRuleId.put(String.valueOf(rule_id), currentElement);
      }//end while
      
      //
      //Read in the dependencies for each rule. Add each dependeny as a "parent" for each rule.
      select = "select rule_id, condition from rule_conditions";
      resultSet = statement.executeQuery(select);
      while(resultSet.next()){
        rule_id = resultSet.getInt("rule_id");
        rule_condition = resultSet.getString("condition");
        parentElement = (ModelElement)elementsByLogicName.get(rule_condition);
        currentElement = (ExtendedModelElement)elementsByRuleId.get(String.valueOf(rule_id));
        parentElement.addChildElement(currentElement);
      }

      //
      //OK, so now we have s set of nodes, indexed by id etc, all connected to
      //each other via parent relationships. We have to attach the ones with no
      //parent to the root.
      Iterator elements = elementsByLogicName.values().iterator();
      while(elements.hasNext()){
        currentElement = (ExtendedModelElement)elements.next();
        if(currentElement.getParentElements().size() <= 0){
          rootElement.addChildElement(currentElement);
        }
      }
      
      //
      //Keep a Collection - as an attribute of the root element - which 
      //contains all the nodes.
      if(getLogger().isLoggingMedium()){
        getLogger().logMedium("Adding "+elementsByLogicName.size()+" to model with key "+PathStepsModel.PATH_STEPS_PANEL_ALL_NODES);
      }

      rootElement.addProperty(PathStepsModel.PATH_STEPS_PANEL_ALL_NODES, elementsByLogicName.values());
      
      rootElement.addProperty(
        PathStepsModel.PATH_STEPS_PANEL_GRAPH_LAYOUT_CONFIGURATION,
        getView().getApplication().readGraphLayoutConfiguration()
      );
      
    }catch(SQLException exception){
      if(getLogger().isLoggingMedium()){
        getLogger().logMedium("SQL Problems reading pipline database", exception);
      }
      
      throw new NonFatalAException("SQL Problems reading pipline database", exception);
    }

   
    if(getLogger().isLoggingHigh()){
      getLogger().logHigh("Datamodel after db read: ");
      getLogger().logHigh(rootElement.toString());
    }
    return rootElement;
  }
  
}
