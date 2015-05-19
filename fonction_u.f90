module fonction_u

  use exact
  
contains

  function simpson(w,c,a,b)
    implicit none
    real*8,intent(in)::w,c,a,b
    real*8::simpson
    simpson=((b-a)/6.)*(phi2(w,c,a)+4*phi2(w,c,(a+b)/2.)+phi2(w,c,b))
  end function simpson

  function integrale_infini (w,c,f,h)
    implicit none
    real*8,intent(in)::w,c,f,h
    real*8::a,b,suivant,integrale_infini
    a=f-h
    b=f
    integrale_infini = simpson(w,c,a,b)
    b=b-h
    a=a-h
    suivant=simpson(w,c,a,b)
    do while (suivant>(h/10.))
       integrale_infini = integrale_infini + suivant
       b=b-h
       a=a-h
       suivant=simpson(w,c,a,b)
    end do
    !print*,a
  end function integrale_infini 

  function integrale_positive (w,c,f,h)
    implicit none
    real*8,intent(in)::w,c,f,h
    real*8::a,b,suivant,integrale_positive
    a=0.0D0
    b=h
    integrale_positive = simpson(w,c,a,b)
    b=b+h
    a=a+h
    suivant=simpson(w,c,a,b)
    do while (f>b)
       integrale_positive = integrale_positive + suivant
       b=b+h
       a=a+h
       suivant=simpson(w,c,a,b)
    end do
  end function integrale_positive

  function integrale_phi2 (w,c,f,h)
    implicit none
    real*8,intent(in)::w,c,f,h
    real*8::integrale_phi2
    if (f>0) then
       integrale_phi2=integrale_infini(w,c,0.0D0,h)+integrale_positive (w,c,f,h)
    else
       integrale_phi2=integrale_infini(w,c,f,h)
    end if
  end function integrale_phi2
    
  subroutine exact_u0 (w,c,x,h,U)
    implicit none
    real*8 ,intent(in)::w,c,h
    real*8 ,dimension(:) ,intent(in)::x
    complex*8 ,dimension(:) ,intent(out)::U
    real*8 ,dimension(:),allocatable::phi
    real*8::var
    integer::i,n
    n=size(x,1)
    allocate(phi(n))
    call exact_phi(w,c,x,phi)
    do i=1,size(x,1)
       var=integrale_phi2 (w,c,x(i),h)
       var=((c*x(i))/2.)-(3./4.)*var
       U(i)=cmplx(phi(i)*cos(var),phi(i)*sin(var))
    end do
    U=phi*U
    deallocate(phi)
  end subroutine exact_u0

 subroutine ctoreal(U,a,b)
    implicit none
    complex*8 ,dimension(:) ,intent(in)::U
    real*8 ,dimension(:) ,intent(out)::a,b
    integer::i
    do i=1,size(U,1)
       a(i)=real(U(i))
       b(i)=aimag(U(i))
    end do
  end subroutine ctoreal
  

    

end module fonction_u
