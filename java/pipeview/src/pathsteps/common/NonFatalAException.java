package pathsteps.common;

public class NonFatalAException extends RuntimeException{
  public NonFatalAException(String message){
    super(message);
  }
  
  public NonFatalAException(String message, Throwable originalException){
    super(message, originalException);
  }
}
