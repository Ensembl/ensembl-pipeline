package pathsteps.test;

import pathsteps.common.*;
import pathsteps.gui.*;
import pathsteps.model.*;
import junit.framework.*;

public class TestRunner {

  Application application;
  PathStepsTestView view;
  PathStepsModel model;
  
  public TestRunner(
    Application application,
    PathStepsTestView view,
    PathStepsModel model
  ){
    this.application = application;
    this.view = view;
    this.model = model;
  }

  public Application getApplication(){
    return application;
  }
  
  public PathStepsModel getModel(){
    return model;
  }
  
  public PathStepsTestView getView(){
    return view;
  }
  
  public void run(){
    TestSuite suite = new TestSuite();
    
    //
    //Which test cases are added to the suite should be configurable - via the test.conf file!
    suite.addTest(new CreateNewLayoutConfigurationDialogTest(this));
    
    TestResult result = new TestResult();
    suite.run(result);
    
    System.err.println(result.toString()+result.wasSuccessful()+" failures: "+result.failureCount()+" errors: "+result.errorCount());
    System.err.println("Failures: ");
    java.util.Enumeration failures = result.failures();
    while(failures.hasMoreElements()){
      TestFailure failure = (TestFailure)failures.nextElement();
      System.err.println(failure);
      failure.thrownException().printStackTrace();
    }
    
    System.err.println("Errors: ");
    java.util.Enumeration errors = result.errors();
    while(errors.hasMoreElements()){
      TestFailure failure = (TestFailure)errors.nextElement();
      System.err.println(failure);
      failure.thrownException().printStackTrace();
    }

  }
}
