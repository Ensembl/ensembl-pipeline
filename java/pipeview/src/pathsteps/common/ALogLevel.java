package pathsteps.common;
import java.util.logging.*;

public class ALogLevel extends Level{

  public static final Level LOW = new ALogLevel("LOW", Level.INFO.intValue());
  public static final Level MEDIUM = new ALogLevel("MEDIUM", Level.FINE.intValue());
  public static final Level HIGH = new ALogLevel("HIGH", Level.FINEST.intValue());
  
  public ALogLevel(String name, int value) {
    super(name, value);
  }

}
