package pathsteps.common;
import java.util.*;
import java.text.*;
import java.util.logging.*;

public class JavaLogger implements ALogger{
  private static String lowIndent = "";
  private static String mediumIndent = "\t";
  private static String highIndent = "\t\t";
  private DateFormat format = SimpleDateFormat.getInstance();
  private String prefix ="";

  //
  //LOW -> INFO
  //MEDIUM -> FINE
  //HIGH -> FINEST
  private Logger logger;

  public JavaLogger(){}
  
  public JavaLogger(String prefix){
    this.prefix = prefix;
    logger = Logger.getLogger("pathsteps.common."+prefix);
  }
  
  private Logger getLogger(){
    return logger;
  }
  
  private void flush(){
    for(int i=0; i<getLogger().getHandlers().length; i++){
      getLogger().getHandlers()[i].flush();
    }
  }
  
  public void logLow(String message){
    ALogRecord record = new ALogRecord(ALogLevel.LOW, message, getPrefix());
    getLogger().log(record);
    flush();
  }
  
  public void logLow(String message, Throwable exception){
    ALogRecord record = new ALogRecord(ALogLevel.LOW, message, getPrefix());
    record.setThrown(exception);
    getLogger().log(record);
    flush();
  }
  
  public boolean isLoggingLow(){
    return getLogger().isLoggable(ALogLevel.LOW);
  }
  
  public void logMedium(String message){
    ALogRecord record = new ALogRecord(ALogLevel.MEDIUM, message, getPrefix());
    getLogger().log(record);
    flush();
  }
  
  public void logMedium(String message, Throwable exception){
    ALogRecord record = new ALogRecord(ALogLevel.MEDIUM, message, getPrefix());
    record.setThrown(exception);
    getLogger().log(record);
    flush();
  }
  
  public boolean isLoggingMedium(){
    return getLogger().isLoggable(ALogLevel.MEDIUM);
  }
  
  public void logHigh(String message){
    ALogRecord record = new ALogRecord(ALogLevel.HIGH, message, getPrefix());
    getLogger().log(record);
    flush();
  }
  
  public void logHigh(String message, Throwable exception){
    ALogRecord record = new ALogRecord(ALogLevel.HIGH, message, getPrefix());
    record.setThrown(exception);
    getLogger().log(record);
    flush();
  }
  
  public boolean isLoggingHigh(){
    return getLogger().isLoggable(ALogLevel.HIGH);
  }
 
  private String getPrefix(){
    return prefix;
  }
}