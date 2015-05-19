module erreur
  use exact
  use tableau
  use resolve

contains
  function norme2(V)
    implicit none
    real*8,dimension(:),intent(in)::V
    real*8::norme2
    integer::i,n
    n=size(V,1)
    norme2=0.0D0
    do i=1,n
       norme2=norme2+V(i)**2
    end do
    norme2=dsqrt(norme2)
  end function norme2

  function normeinfinie(V)
    implicit none
    real*8,dimension(:),intent(in)::V
    real*8::normeinfinie
    integer::i,n
    n=size(V,1)
    normeinfinie=abs(V(1))
    do i=2,n
       if (normeinfinie<abs(V(i))) then
          normeinfinie=abs(V(i))
       end if
    end do
  end function normeinfinie

  subroutine erreur_Euler_e(a,b,h,w,c,n2,nf)
    implicit none
    real*8,intent(in)::a,b,h,w,c
    real*8,intent(out)::n2,nf
    real*8,dimension(:),allocatable::x,y,ye,e
    real*8,dimension(2)::M,M_next
    integer::n,i,k
    call cree_tableau_pas(x,a,b,h)
    n=size(x,1)
    allocate(y(n))
    allocate(ye(n))
    call exact_phi(w,c,x,y)   
    M(1)=y(1)
    M(2)=exact_Dphi(w,c,x(1))   
    ye(1)=M(1)
    do i=2,n
       call Euler_e(2,M,M_next,h,w,c)
       ye(i)=M_next(1)
       M=M_next
    end do
    allocate(e(n))
    e=y-ye
    n2=norme2(e)/norme2(y)
    nf=normeinfinie(e)/normeinfinie(y)
    deallocate(x)
    deallocate(y)
    deallocate(ye)
    deallocate(e)
  end subroutine erreur_Euler_e


 subroutine erreur_Euler_i(a,b,h,w,c,n2,nf)
    implicit none
    real*8,intent(in)::a,b,h,w,c
    real*8,intent(out)::n2,nf
    real*8,dimension(:),allocatable::x,y,ye,e
    real*8,dimension(2)::M,M_next
    integer::n,i,k
    call cree_tableau_pas(x,a,b,h)
    n=size(x,1)
    allocate(y(n))
    allocate(ye(n))
    call exact_phi(w,c,x,y)   
    M(1)=y(1)
    M(2)=exact_Dphi(w,c,x(1))   
    ye(1)=M(1)
    do i=2,n
       call Euler_i(2,M,M_next,h,w,c)
       ye(i)=M_next(1)
       M=M_next
    end do
    allocate(e(n))
    e=y-ye
    n2=norme2(e)/norme2(y)
    nf=normeinfinie(e)/normeinfinie(y)
    deallocate(x)
    deallocate(y)
    deallocate(ye)
    deallocate(e)
  end subroutine erreur_Euler_i
  
  subroutine erreur_RK4(a,b,h,w,c,n2,nf)
    implicit none
    real*8,intent(in)::a,b,h,w,c
    real*8,intent(out)::n2,nf
    real*8,dimension(:),allocatable::x,y,yk,e
    real*8,dimension(2)::M,M_next
    integer::n,i,k
    call cree_tableau_pas(x,a,b,h)
    n=size(x,1)
    allocate(y(n))
    allocate(yk(n))
    call exact_phi(w,c,x,y)   
    M(1)=y(1)
    M(2)=exact_Dphi(w,c,x(1))   
    yk(1)=M(1)
    do i=2,n
       call RK4(2,M,M_next,h,w,c)
       yk(i)=M_next(1)
       M=M_next
    end do
    allocate(e(n))
    e=y-yk
    n2=norme2(e)/norme2(y)
    nf=normeinfinie(e)/normeinfinie(y)
    deallocate(x)
    deallocate(y)
    deallocate(yk)
    deallocate(e)    
  end subroutine erreur_RK4
  
  subroutine erreur_Adam4(a,b,h,w,c,n2,nf)
    implicit none
    real*8,intent(in)::a,b,h,w,c
    real*8,intent(out)::n2,nf
    real*8,dimension(:),allocatable::x,y,yk,e
    real*8,dimension(2)::M1,M2,M3,M4,M_next
    integer::n,i,k
    call cree_tableau_pas(x,a,b,h)
    n=size(x,1)
    allocate(y(n))
    allocate(yk(n))
    call exact_phi(w,c,x,y)   
    M4(1)=y(1)
    M4(2)=exact_Dphi(w,c,x(1)) 
    M3(1)=y(2)
    M3(2)=exact_Dphi(w,c,x(2)) 
    M2(1)=y(3)
    M2(2)=exact_Dphi(w,c,x(3)) 
    M1(1)=y(4)
    M1(2)=exact_Dphi(w,c,x(4)) 
    yk(1)=M4(1)
    yk(2)=M3(1)
    yk(3)=M2(1)
    yk(4)=M1(1)
    do i=5,n
       call Adam4(2,M1,M2,M3,M4,M_next,h,w,c)
       yk(i)=M_next(1)
       M4=M3
       M3=M2    
       M2=M1
       M1=M_next  
    end do
    allocate(e(n))
    e=y-yk
    n2=norme2(e)/norme2(y)
    nf=normeinfinie(e)/normeinfinie(y)
    deallocate(x)
    deallocate(y)
    deallocate(yk)
    deallocate(e)    
  end subroutine erreur_Adam4
  
  subroutine erreur_Adam4cor(a,b,h,w,c,n2,nf)
    implicit none
    real*8,intent(in)::a,b,h,w,c
    real*8,intent(out)::n2,nf
    real*8,dimension(:),allocatable::x,y,yk,e
    real*8,dimension(2)::M1,M2,M3,M4,M_next
    integer::n,i,k
    call cree_tableau_pas(x,a,b,h)
    n=size(x,1)
    allocate(y(n))
    allocate(yk(n))
    call exact_phi(w,c,x,y)   
    M4(1)=y(1)
    M4(2)=exact_Dphi(w,c,x(1)) 
    M3(1)=y(2)
    M3(2)=exact_Dphi(w,c,x(2)) 
    M2(1)=y(3)
    M2(2)=exact_Dphi(w,c,x(3)) 
    M1(1)=y(4)
    M1(2)=exact_Dphi(w,c,x(4)) 
    yk(1)=M4(1)
    yk(2)=M3(1)
    yk(3)=M2(1)
    yk(4)=M1(1)
    do i=5,n
       call Adam4cor(2,M1,M2,M3,M4,M_next,h,w,c)
       yk(i)=M_next(1)
       M4=M3
       M3=M2    
       M2=M1
       M1=M_next  
    end do
    allocate(e(n))
    e=y-yk
    n2=norme2(e)/norme2(y)
    nf=normeinfinie(e)/normeinfinie(y)
    deallocate(x)
    deallocate(y)
    deallocate(yk)
    deallocate(e)    
  end subroutine erreur_Adam4cor

end module erreur
