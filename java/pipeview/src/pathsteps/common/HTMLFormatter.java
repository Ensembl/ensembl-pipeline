package pathsteps.common;

import java.util.logging.*;
import java.util.*;
import java.text.*;

public class HTMLFormatter extends Formatter{

  public static String INDENT="................";
//  public static String INDENT="";
  public static String SIZE_PREFIX ="<small>";
  public static String SIZE_SUFFIX ="</small>";
  
  public HTMLFormatter() {
  }

  public String format(LogRecord logRecord) {
    ALogRecord record =  null;
    String indent = "";
    String sizePrefix = "";
    String sizeSuffix = "";
    StringBuffer buffer = new StringBuffer(1000);
    
    if(!(logRecord instanceof ALogRecord)){
      throw new FatalAException(
        "pathsteps.common.HTMLFormatter can only be used to format instances of pathsteps.common.ALogRecord"
      );
    }
    
    record = (ALogRecord)logRecord;
    
    if(record.getDomain().equals(ALogRecord.APP)){
      indent = INDENT;
    }else if(record.getDomain().equals(ALogRecord.CONTROL)){
      indent += INDENT+INDENT;
    }else if(record.getDomain().equals(ALogRecord.MODEL)){
      indent += INDENT+INDENT+INDENT;
    }
    
    if (record.getLevel() == ALogLevel.MEDIUM) {
      sizePrefix += SIZE_PREFIX;
      sizeSuffix = SIZE_SUFFIX+sizeSuffix;
    }else if (record.getLevel() == ALogLevel.HIGH) {
      sizePrefix += SIZE_PREFIX+SIZE_PREFIX;
      sizeSuffix = SIZE_SUFFIX+SIZE_SUFFIX+sizeSuffix;
    }
    
    buffer
      .append(indent)
      .append(sizePrefix)
      .append(formatMessage(logRecord))
      .append(sizeSuffix)
      .append("<br>")
      .append('\n');
    
    return buffer.toString();
  }
  
  // This method is called just after the handler using this
  // formatter is created
  public String getHead(Handler handler) {
    return "<HTML><HEAD>"+(new Date())+"</HEAD><BODY><BASEFONT FACE=\"ARIAL, VERDANA, COURIER\", SIZE=4>\n";
  }

  // This method is called just after the handler using this
  // formatter is closed
  public String getTail(Handler handler) {
    return "</BODY></HTML>\n";
  }    
}
