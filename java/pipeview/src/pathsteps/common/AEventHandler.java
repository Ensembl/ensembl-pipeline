package pathsteps.common;

public interface AEventHandler{
  public Application getApplication();
  public void setApplication(Application application);
  public ALogger getLogger();
  public void setLogger(ALogger logger);
  public void initialise();
  public void notifyEventForKey(String key);
  public void doActionForKey(String key);
}
