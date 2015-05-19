module resolveLUComplex


contains
  subroutine prodvect(M,N,P)
    implicit none
    complex*8,dimension(:),intent(in)::M,N
    complex*8,dimension(:,:),intent(out)::P
    integer::i,j,h,f,k
    h=size(M,1)
    f=size(N,1)
    do i=1,h
       do j=1,f
          P(i,j)=M(i)*N(j)
       end do
    end do
  end subroutine prodvect


  subroutine LU_decomposition(A,L,U)
    implicit none
    complex*8,dimension(:,:),intent(in)::A
    complex*8,dimension(:,:),intent(out)::L,U
    complex*8,dimension(:,:),allocatable::M,P
    complex*8,dimension(:),allocatable::P1,P2
    integer::i,n
    n=size(A,1)
    allocate(M(n,n))
    M=A
    do i=1,n
       M(i,i+1:n)=M(i,i+1:n)/M(i,i)
       allocate(P1(n-(i+1)))
       allocate(P2(n-(i+1)))
       allocate(P(n-(i),n-(i)))
       P1=M(i+1:n,i)
       P2=M(i,i+1:n)
       call prodvect(P1,P2,P)
       M(i+1:n,i+1:n)=M(i+1:n,i+1:n)-P!-M(i+1:n,i)*M(i,i+1:n)
       deallocate(P)
       deallocate(P1)
       deallocate(P2)
    end do
    L=cmplx(0.0D0,0.0D0)
    U=cmplx(0.0D0,0.0D0)
    do i=1,n
       L(1:i,i)=M(1:i,i)
       U(i:n,i)=M(i:n,i)
       L(i,i)=cmplx(1.0D0,0.0D0)
    end do
    deallocate(M)
  end subroutine LU_decomposition
  

  subroutine resolve_sysLU (L,U,b,x)
    implicit none
    complex*8,dimension(:,:),intent(in)::L,U
    complex*8,dimension(:),intent(in)::b
    complex*8,dimension(:),intent(out)::x
    complex*8,dimension(:),allocatable::y
    complex*8::var
    integer::i,j,n
    n=size(L,1)
    allocate(y(n))
    !resolve Ly=b
    y(1)=b(1)    
    do i=2,n
       var=cmplx(0.0D0,0.0D0)
       do j=1,i-1
          var=var+L(j,i)*y(j)
       end do
       y(i)=b(i)-var
    end do
    !resolve Ux=y
    x(n)=y(n)/U(n,n)
    do i=n-1,1,-1
       var=cmplx(0.0D0,0.0D0)
       do j=i+1,n
          var=var+U(j,i)*x(j)
       end do
       x(i)=(y(i)-var)/U(i,i)
    end do
    deallocate(y)
  end subroutine resolve_sysLU

  subroutine test()
    implicit none 
    complex*8,dimension(4,4)::A,L,U
    complex*8,dimension(4)::b,x
    A(:,1)=(/4.,-1.,-1.,0./)
    A(:,2)=(/-1.,4.,-1.,-1./)
    A(:,3)=(/-1.,-1.,4.,-1./)
    A(:,4)=(/0.,-1.,-1.,4./)
    B=(/2.,1.,1.,2./)
    print*,"-----------------------TEST-----------------------------"
    print*,"A="
    print*,A(:,1)
    print*,A(:,2)
    print*,A(:,3)
    print*,A(:,4)
    print*,"b=",b
    call LU_decomposition(A,L,U)
    call resolve_sysLU (L,U,b,x)
    print*,"L="
    print*,L(:,1)
    print*,L(:,2)
    print*,L(:,3)
    print*,L(:,4)
    print*,"U="
    print*,U(:,1)
    print*,U(:,2)
    print*,U(:,3)
    print*,U(:,4)
    print*,"x=",x
  end subroutine test
end module resolveLUComplex
