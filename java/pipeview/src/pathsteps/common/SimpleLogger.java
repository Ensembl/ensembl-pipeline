package pathsteps.common;
import java.util.*;
import java.text.*;

public class SimpleLogger implements ALogger{
  private static String lowIndent = "";
  private static String mediumIndent = "\t";
  private static String highIndent = "\t\t";
  private DateFormat format = SimpleDateFormat.getInstance();
  private String prefix;

  public SimpleLogger(){}
  
  public SimpleLogger(String prefix){
    this.prefix = prefix;
  }
  
  public void logLow(String message){
    System.out.println(createLowMessage(message));
  }
  
  public void logLow(String message, Throwable exception){
    System.out.println(createLowMessage(message));
    exception.printStackTrace();
  }
  
  public boolean isLoggingLow(){
    return true;
  }
  
  public void logMedium(String message){
    System.out.println(createMediumMessage(message));
  }
  
  public void logMedium(String message, Throwable exception){
    System.out.println(createMediumMessage(message));
    exception.printStackTrace();
  }
  
  public boolean isLoggingMedium(){
    return true;
  }
  
  public void logHigh(String message){
    System.out.println(createHighMessage(message));
  }
  
  public void logHigh(String message, Throwable exception){
    System.out.println(createHighMessage(message));
    exception.printStackTrace();
  }
  
  public boolean isLoggingHigh(){
    return true;
  }
  
  private String getLowIndent(){
    return lowIndent;
  }
  
  private String getMediumIndent(){
    return mediumIndent;
  }
  
  private String getHighIndent(){
    return highIndent;
  }
  
  private String getApplicationTime(){
    return getFormat().format(new Date(System.currentTimeMillis()));
  }
  
  private DateFormat getFormat(){
    return format;
  }
  
  private String getPrefix(){
    return prefix;
  }
  
  private String createLowMessage(String message){
    return getPrefix()+":"+getApplicationTime()+":"+message;
  }
  
  private String createMediumMessage(String message){
    return getMediumIndent()+getPrefix()+":"+getApplicationTime()+":"+message;
  }
  
  private String createHighMessage(String message){
    return getHighIndent()+getPrefix()+":"+getApplicationTime()+":"+message;
  }
}
