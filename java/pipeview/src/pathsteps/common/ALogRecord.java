package pathsteps.common;

import java.util.logging.*;

public class ALogRecord extends LogRecord {
  
  public static String VIEW = "VIEW";
  public static String APP = "APP";
  public static String CONTROL = "CONTROL";
  public static String MODEL = "MODEL";
  
  String domain;
  
  public ALogRecord(Level level, String message, String domain) {
    super(level, message);
    this.domain = domain;
  }
  
  public String getDomain(){
    return domain;
  }
}
