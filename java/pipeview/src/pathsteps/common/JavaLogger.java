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
  
  public void logLow(String message){
    getLogger().info(message);
  }
  
  public void logLow(String message, Throwable exception){
    logLow(message);
    getLogger().throwing(getClass().getName(), "logLow", exception);
  }
  
  public boolean isLoggingLow(){
    if(getLogger().getLevel() != null){
      return getLogger().getLevel().equals(Level.INFO);
    }else if(getLogger().getParent() != null){
      return getLogger().getParent().getLevel().equals(Level.INFO);
    }else{
      return false;
    }
  }
  
  public void logMedium(String message){
    getLogger().fine(message);
  }
  
  public void logMedium(String message, Throwable exception){
    logMedium(message);
    getLogger().throwing(getClass().getName(), "logMedium", exception);
  }
  
  public boolean isLoggingMedium(){
    if(getLogger().getLevel() != null){
      return getLogger().getLevel().equals(Level.FINE);
    }else if(getLogger().getParent() != null){
      return getLogger().getParent().getLevel().equals(Level.FINE);
    }else{
      return false;
    }
  }
  
  public void logHigh(String message){
    getLogger().finest(message);
  }
  
  public void logHigh(String message, Throwable exception){
    logHigh(message);
    getLogger().throwing(getClass().getName(), "logHigh", exception);
  }
  
  public boolean isLoggingHigh(){
    if(getLogger().getLevel() != null){
      return getLogger().getLevel().equals(Level.FINEST);
    }else if(getLogger().getParent() != null){
      return getLogger().getParent().getLevel().equals(Level.FINEST);
    }else{
      return false;
    }
  }
}
