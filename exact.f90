module exact

contains
  subroutine exact_phi(w,c,x,phi)
    implicit none
    real*8 ,intent(in)::w,c
    real*8 ,dimension(:) ,intent(in)::x
    real*8 ,dimension(:) ,intent(out)::phi
    real*8 ::A,B,D
    A=4*w-c**2
    B=dsqrt(w)/(A)
    D=c/(2*dsqrt(w))
    phi=1./dsqrt(B*(dcosh(dsqrt(A)*x)-D))
  end subroutine exact_phi

  function exact_Dphi(w,c,x)
    implicit none
    real*8 ,intent(in)::w,c
    real*8 ,intent(in)::x
    real*8 ::exact_Dphi
    real*8 ::A,B,D
    A=4*w-c**2
    B=dsqrt(w)/(A)
    D=c/(2*dsqrt(w))
    exact_Dphi=(-1./2.)*(B*dsqrt(A)*dsinh(dsqrt(A)*x))*(B*(dcosh(dsqrt(A)*x)-D))**(-3./2.)
  end function exact_Dphi

  function phi2(w,c,x)
    implicit none
    real*8 ,intent(in)::w,c
    real*8 ,intent(in)::x
    real*8::phi2
    real*8 ,dimension(1)::phi,xtab
    xtab(1)=x
    call exact_phi(w,c,xtab,phi)
    phi2=phi(1)**2
  end function phi2

end module exact
