module equation
  
contains

  subroutine fonct(X, FX,w,c)
    implicit none
    real*8, intent(in), dimension (2) :: X
    real*8, intent(out), dimension (2) :: FX
    real*8,intent(in) ::w,c
    real*8  :: a,b,e
    FX=0.0D0
    
    FX(1)=X(2)
    a=(w-((c**2)/4.))*X(1)
    b=(c/2.)*((dabs(X(1)))**2)*X(1)
    e=(-3./16.)*((dabs(X(1)))**4)*X(1)
    FX(2)=a+b+e
  end subroutine fonct
    
end module equation
