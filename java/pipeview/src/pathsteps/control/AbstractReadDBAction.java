package pathsteps.control;

import pathsteps.common.*;
import pathsteps.model.*;

import java.util.*;
import java.sql.*;

/**
 * Common function to whether you're reading the db for the first time, or refreshing
 * information.
**/
public abstract class AbstractReadDBAction extends AAction{

  public AbstractReadDBAction(AEventHandler eventHandler){
    super(eventHandler);
  }
  
  protected void doAction(AView view, AModel aModel){
    PathStepsModel model = (PathStepsModel)aModel;
    Collection databases = null;
    String selectedDatabase = null;
    ModelElement dialogModel = model.getRootElement().getChildElement(model.READ_DB_DIALOG);
    ModelElement layoutDialogModel = model.getRootElement().getChildElement(model.LAYOUT_DIALOG);
    String showDetail = null;
    
    view.getApplication().readHistory(); //picks up the user's last choices from history instead of internal state.
    showDetail = getStringFromHistoryOrConfig(model.LAYOUT_DIALOG_SHOW_JOB_DETAIL);
    if(showDetail == null){
      showDetail = Boolean.FALSE.toString();
    }
    
    Collection list;
    String host = null;
    String port = null;
    String user = null;
    String password = null;
    
    if(isLoggingLow()){
      logLow("Starting ReadSpecifiedDBAction");
    }
    
    if(dialogModel == null){
      throw new FatalAException("ReadDBDialog expecting a model called: "+model.READ_DB_DIALOG+" but none was found");
    }
    
    host = (String)dialogModel.getProperty(model.READ_DB_DIALOG_HOST);
    port = (String)dialogModel.getProperty(model.READ_DB_DIALOG_PORT);
    user = (String)dialogModel.getProperty(model.READ_DB_DIALOG_USER);

    selectedDatabase = (String)dialogModel.getProperty(model.READ_DB_DIALOG_PIPELINE_DB);
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
      createPathStepsPanelModel(selectedDatabase, showDetail)
    );
    
    //System.out.println(model.toString());
    Properties history = view.getApplication().getHistory();
    history.setProperty(model.READ_DB_DIALOG_HOST, host);
    history.setProperty(model.READ_DB_DIALOG_PORT, port);
    history.setProperty(model.READ_DB_DIALOG_USER, user);
    history.setProperty(model.READ_DB_DIALOG_PIPELINE_DB, selectedDatabase);
    
    getView().getApplication().writeHistory(history);
    
    if(isLoggingLow()){
      logLow("Wrote application history");
    }
  }

  protected ModelElement createPathStepsPanelModel(String selectedDatabase, String showJobDetail){
    
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
    
    if(isLoggingMedium()){
      logMedium("Started creating model");
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
        currentElement.addProperty(ModelElement.SHOW_DETAIL, showJobDetail);
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
        if(currentElement != null){
          parentElement.addChildElement(currentElement);
        }
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
      if(isLoggingMedium()){
        logMedium("Adding "+elementsByLogicName.size()+" to model with key "+PathStepsModel.PATH_STEPS_PANEL_ALL_NODES);
      }

      addMonitorInformation(statement, elementsByLogicName.values());

      rootElement.addProperty(PathStepsModel.PATH_STEPS_PANEL_ALL_NODES_MAP, elementsByLogicName);
      rootElement.addProperty(PathStepsModel.PATH_STEPS_PANEL_ALL_NODES, elementsByLogicName.values());

      if(isLoggingMedium()){
        logMedium("Adding graph layout configuration properties object to model");
      }

      rootElement.addProperty(
        PathStepsModel.PATH_STEPS_PANEL_GRAPH_LAYOUT_CONFIGURATION,
        getView().getApplication().readGraphLayoutConfiguration()
      );
      
      rootElement.addProperty(
        PathStepsModel.PATH_STEPS_PANEL_SHOW_JOB_DETAIL,
        showJobDetail
      );

      statement.close();
      conn.close();
      
    }catch(SQLException exception){
      if(isLoggingMedium()){
        logMedium("SQL Problems reading pipline database", exception);
      }
      
      throw new NonFatalAException("SQL Problems reading pipline database", exception);
    }
   
    if(isLoggingHigh()){
      logHigh("Datamodel after db read: ");
      logHigh(rootElement.toString());
    }
    return rootElement;
  }

  /**
   * Reads the pipeline tables, adds job-completion etc information as attributes
   * onto each ModelElement passed in via the collection
  **/
  protected void addMonitorInformation(
    Statement statement, 
    Collection allElements
  ) throws SQLException{

    Iterator elements = allElements.iterator();
    HashMap currentlyRunningStatuses = getCurrentlyRunningJobStatuses(statement);
    HashMap finishedCounts = getFinishedCounts(statement);
    ModelElement element;
    String elementName;
    Map runningStatusMap;

    while(elements.hasNext()){
      element = (ModelElement)elements.next();
      elementName = element.getKey();
      if(finishedCounts.get(elementName)!= null){
        element.addProperty(ModelElement.FINISHED_COUNT, finishedCounts.get(elementName));
      }
      
      if(currentlyRunningStatuses.get(elementName)!=null){
        runningStatusMap = (Map)currentlyRunningStatuses.get(elementName);

        if(runningStatusMap.get(ModelElement.SUBMITTED) != null){
          element.addProperty(
            ModelElement.SUBMITTED, 
            runningStatusMap.get(ModelElement.SUBMITTED)
          );
        }
        
        if(runningStatusMap.get(ModelElement.RUNNING) != null){
          element.addProperty(
            ModelElement.RUNNING, 
            runningStatusMap.get(ModelElement.RUNNING)
          );
        }
        
        if(runningStatusMap.get(ModelElement.FAILED) != null){
          element.addProperty(
            ModelElement.FAILED, 
            runningStatusMap.get(ModelElement.FAILED)
          );
        }
      }
    }

  }

  /**
   * Return a Hash-of-Hashes. The first hash is keyed by analysis name. Each value
   * in the first hash is another hash, which is keyed by job status. Each value in
   * the second hash is an (Integer) count of the number of jobs in that status.
   * This counts all unfinished jobs by analysis id.
  **/
  protected HashMap getCurrentlyRunningJobStatuses(Statement statement) throws SQLException{
    String analysis_name;
    String old_analysis_name;
    String job_status;
    int count;
    HashMap returnMap = new HashMap();
    HashMap countsForAnalysis = null;
    String select = 
      (new StringBuffer())
        .append("select analysis.logic_name, job_status.status, count(*)  ")
        .append("from analysis, job_status, job ")
        .append("where job.job_id = job_status.job_id and ")
        .append("analysis.analysis_id = job.analysis_id and  ")
        .append("job_status.is_current = 'y' ")
        .append("group by analysis.logic_name, job_status.status ")
        .append("order by analysis.logic_name")
        .toString(); 

    ResultSet resultSet = statement.executeQuery(select);
    old_analysis_name = "";
    
    while(resultSet.next()){
      analysis_name = resultSet.getString(1);
      job_status = resultSet.getString(2);
      count = resultSet.getInt(3);

      if(
        job_status.equals(ModelElement.SUBMITTED) ||
        job_status.equals(ModelElement.FAILED) ||
        job_status.equals(ModelElement.RUNNING)
      ){
        if(!old_analysis_name.equals(analysis_name)){
          countsForAnalysis = new HashMap();
          returnMap.put(analysis_name, countsForAnalysis);
        }
        if(job_status.equals(ModelElement.RUNNING)){
          countsForAnalysis.put(ModelElement.RUNNING, new Integer(count));
        }else if(job_status.equals(ModelElement.FAILED)){
          countsForAnalysis.put(ModelElement.FAILED, new Integer(count));
        }else if(job_status.equals(ModelElement.SUBMITTED)){
          countsForAnalysis.put(ModelElement.SUBMITTED, new Integer(count));
        }
      }
      
      old_analysis_name = analysis_name;
    }
    return returnMap;
  }
  
  /**
   * Return a hashmap keyed by analysis name, of the count of finished jobs for
   * each analysis id.
  **/
  private HashMap getFinishedCounts(Statement statement) throws SQLException{  
    HashMap returnMap = new HashMap();
    String analysis_name = null;
    int count;
    String select = 
      (new StringBuffer())
        .append("select analysis.logic_name, count(*) ")
        .append("from input_id_analysis, analysis ")
        .append("where analysis.analysis_id = input_id_analysis.analysis_id ")
        .append(" group by analysis.analysis_id").toString();

    ResultSet resultSet = statement.executeQuery(select);
    while(resultSet.next()){
      analysis_name = resultSet.getString(1);
      count = resultSet.getInt(2);
      returnMap.put(analysis_name, new Integer(count));
    }

    return returnMap;
  }
  
/*
  private ModelElement createFakePathStepsModel(String selectedDatabase){
    
    ModelElement rootElement = new ModelElement(PathStepsModel.PATH_STEPS_PANEL);

    //ModelElement submitContig = rootElement.createChildElement("SubmitContig");
    //ModelElement repeatMask = submitContig.createChildElement("RepeatMask");
    //ModelElement CpG = submitContig.createChildElement("CpG");
    //ExtendedModelElement accumulator = new ExtendedModelElement("Accumulator");
    //repeatMask.addChildElement(accumulator);
    //CpG.addChildElement(accumulator);

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
*/    
}
