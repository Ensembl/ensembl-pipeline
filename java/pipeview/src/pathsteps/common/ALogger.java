package pathsteps.common;

public interface ALogger{
  public void logLow(String message);
  public void logLow(String message, Throwable exception);
  public boolean isLoggingLow();
  public void logMedium(String message);
  public void logMedium(String message, Throwable exception);
  public boolean isLoggingMedium();
  public void logHigh(String message);
  public void logHigh(String message, Throwable exception);
  public boolean isLoggingHigh();
}
