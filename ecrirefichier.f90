module ecrirefichier
  
contains
  subroutine ecrirefichier2d(name,sta,x,y)
    implicit none
    real*8 ,dimension(:) ,intent(in)::x,y
    character(len=*),intent(in)::name,sta
    integer::i,n
    n=size(x,1)
    open(unit=10,file=name,status=sta)
    do i=1,n
       write(10,*) x(i), y(i)
    end do
    close(10)
    print*,name,' crée'
  end subroutine ecrirefichier2d

  subroutine ecrirefichier3d(name,sta,x,y,z)
    implicit none
    real*8 ,dimension(:) ,intent(in)::x,y,z
    character(len=*),intent(in)::name,sta
    integer::i,n
    n=size(x,1)
    open(unit=10,file=name,status=sta)
    do i=1,n
       write(10,*) x(i), y(i), z(i)
    end do
    close(10)
    print*,name,' crée'
  end subroutine ecrirefichier3d

end module ecrirefichier
