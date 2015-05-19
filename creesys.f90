module systeme


contains
  
  subroutine cree_matriceA (A,U,dt,dx)
    implicit none
    complex*8,dimension(:,:),intent(out)::A
    complex*8,dimension(:),intent(in)::U
    real*8,intent(in)::dt,dx
    integer::i,n
    complex*8::j
    j=cmplx(0,1.)
    n=size(U,1)
    A=cmplx(0.0D0,0.0D0)
    A(1,1)=(j/dt)-(1./(dx**2))
    A(1,2)=(1./(2.*dx**2))+((j*(abs(U(2))**2))/(4.*dx))
    do i=2,n-1
       A(i,1+i)=(1./(2.*dx**2))+((j*(abs(U(i+1))**2))/(4.*dx))
       A(i,i)=(j/dt)-(1./(dx**2))
       A(i,i-1)=(1./(2.*dx**2))+((j*(abs(U(i-1))**2))/(4.*dx))
    end do
    A(n,n)=(j/dt)-(1./(dx**2))
    A(n,n-1)=(1./(2.*dx**2))+((j*(abs(U(n-1))**2))/(4.*dx))
  end subroutine cree_matriceA

  subroutine cree_vecteurB (B,U,dt,dx)
    implicit none
    complex*8,dimension(:),intent(out)::B
    complex*8,dimension(:),intent(in)::U
    real*8,intent(in)::dt,dx
    integer::i,n
    complex*8::j,b1,b2,b3
    j=cmplx(0,1)
    n=size(U,1)
    B=cmplx(0.0D0,0.0D0)
    B(1)=((j/dt)+(1./(dx**2)))*U(1)-(((j*(abs(U(2))**2))/(4.*dx))+(1./(2.*dx**2)))*U(2)
    do i=2,n-1
       b1=(((j*(abs(U(i-1))**2))/(4.*dx)))-(1./(2.*dx**2))
       b2=((j/dt)+(1./(dx**2)))
       b3=-(((j*(abs(U(i+1))**2))/(4.*dx)))-(1./(2.*dx**2))
       B(i)=b1*U(i-1)+b2*U(i)+b3*U(i+1)
    end do
       b1=(((j*(abs(U(n-1))**2))/(4.*dx)))-(1./(2.*dx**2))*U(n-1)
       b2=((j/dt)+(1./(dx**2)))*U(n)
       B(n)=b1+b2
  end subroutine cree_vecteurB




end module systeme
