package pathsteps.common;

public interface AEventHandler{
  public Application getApplication();
  public void setApplication(Application application);
  public ALogger getLogger();
  public void setLogger(ALogger logger);
  public void initialise();
  public void notifyEventForKey(String key);
  public void doActionForKey(String key);
    
  public boolean isLoggingLow();
  public void logLow(String message);
  public void logLow(String message, Throwable exception);
  
  public boolean isLoggingMedium();
  public void logMedium(String message);
  public void logMedium(String message, Throwable exception);
  
  public boolean isLoggingHigh();
  public void logHigh(String message);
  public void logHigh(String message, Throwable exception);  
}
