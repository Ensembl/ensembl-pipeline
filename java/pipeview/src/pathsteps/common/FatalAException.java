package pathsteps.common;

public class FatalAException extends RuntimeException{
  public FatalAException(String message){
    super(message);
  }
  
  public FatalAException(String message, Throwable originalException){
    super(message, originalException);
  }
}
