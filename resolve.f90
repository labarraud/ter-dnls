module resolve
use equation

contains
 subroutine Euler_e(m, X, X_next, dt,w,c)
    
   implicit none
   
   integer, intent(in) :: m
   real*8, intent(in), dimension(m) :: X
   real*8, intent(out), dimension(m) :: X_next
   real*8, dimension(m) :: FX
    real*8,intent(in)::w,c,dt
   call fonct(X,FX,w,c)
   X_next=X+dt*FX
    
 end subroutine Euler_e

 subroutine Euler_i(m, X, X_next, dt,w,c)
   implicit none
   integer, intent(in) :: m
   real*8, intent(in), dimension(m) :: X
   real*8, intent(out), dimension(m) :: X_next
   real*8, dimension(m) :: FX
   real*8,intent(in)::w,c,dt
   call Euler_e(m, X, X_next, dt,w,c)
   call fonct(X_next,FX,w,c)
   X_next=X+dt*FX
 end subroutine Euler_i

 subroutine RK4(m, X, X_next, dt,w,c)
    implicit none
    integer, intent(in) :: m
    real*8,intent(in)::w,c,dt
    real*8, intent(in), dimension(m) :: X
    real*8, intent(out), dimension(m) :: X_next
    real*8, dimension(m) :: FX
    real*8,dimension(m) :: k1,k2,k3,k4
    call fonct(X,FX,w,c)
    k1=FX
    call fonct(X+(dt/2)*k1,FX,w,c)
    k2=FX
    call fonct(X+(dt/2)*k2,FX,w,c)
    k3=FX
    call fonct(X+dt*k3,FX,w,c)
    k4=FX
    X_next=X+(dt/6.)*(k1+2*k2+2*k3+k4)
  end subroutine RK4
 
 subroutine Adam4(m, X1,X2,X3,X4, X_next, dx,w,c)
    implicit none
    integer, intent(in) :: m
    real*8,intent(in)::w,c,dx
    real*8, intent(in), dimension(m) :: X1,X2,X3,X4
    real*8, intent(out), dimension(m) :: X_next
    real*8, dimension(m) :: FX1,FX2,FX3,FX4
    call fonct(X1,FX1,w,c)
    call fonct(X2,FX2,w,c)
    call fonct(X3,FX3,w,c)
    call fonct(X4,FX4,w,c)
    X_next=X1+(dx/24.)*(55*FX1-59*FX2+37*FX3-9*FX4)
  end subroutine Adam4


!methode adam mouthon
 subroutine Adam4cor(m, X1,X2,X3,X4, X_next, dx,w,c)
    implicit none
    integer, intent(in) :: m
    real*8,intent(in)::w,c,dx
    real*8, intent(in), dimension(m) :: X1,X2,X3,X4
    real*8, intent(out), dimension(m) :: X_next
    real*8, dimension(m) :: FX1,FX2,FX3,FX,X
    call Adam4(m, X1,X2,X3,X4,X, dx,w,c)
    call fonct(X1,FX1,w,c)
    call fonct(X2,FX2,w,c)
    call fonct(X3,FX3,w,c)
    call fonct(X,FX,w,c)
    X_next=X1+(dx/24.)*(FX3-5*FX2+19*FX1+9*FX)
  end subroutine Adam4cor


end module resolve
