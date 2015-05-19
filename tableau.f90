module tableau
  
contains
  subroutine cree_tableau_pas(T,a,b,h)
    implicit none
    real*8 ,intent(in)::a,b,h
    real*8 ,dimension(:) ,allocatable,intent(inout)::T
    integer::n,i
    real*8::z
    !n=floor((b-a)/h)+1
    n=1
    z=a
    do while (z<b)
       z=z+h
       n=n+1
    end do
    allocate(T(n))
    T(1)=a
    do i=2,n 
       T(i)=T(i-1)+h
    end do
  end subroutine cree_tableau_pas

subroutine cree_tableau_pas_trou(T,a,b,c,d,h)
    implicit none
    real*8 ,intent(in)::a,b,c,d,h
    real*8 ,dimension(:) ,allocatable,intent(inout)::T
    integer::n,i,m
    real*8::z
    !n=floor((b-a)/h)+1
    n=1
    z=a
    do while (z<b)
       z=z+h
       n=n+1
    end do
    m=n
    n=n+1
    z=c
    do while (z<d)
       z=z+h
       n=n+1
    end do
    allocate(T(n))
    T(1)=a
    do i=2,m
       T(i)=T(i-1)+h
    end do
    T(m+1)=c
    do i=m+2,n
       T(i)=T(i-1)+h
    end do
  end subroutine cree_tableau_pas_trou

    
end module tableau
