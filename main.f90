program main   

  use exact
  use ecrirefichier
  use tableau
  use resolve
  use erreur
  use fonction_u
  use systeme
  use resolveLUComplex

  implicit none   
  real*8,dimension(:),allocatable::x,y,ye,yei,yk,e_n2,RK4_n2,e_nf,RK4_nf,htab,adam4_n2,adam4_nf,ya,ei_n2,ei_nf,ac_n2,ac_nf,r,im
  complex*8,dimension(:,:),allocatable::A1,L,U
  complex*8,dimension(:),allocatable::U0,B1,U1
  real*8,dimension(2)::M,M_next,M1,M2,M3,M4  
  !real*8,dimension(5)::htab 
  real*8::a,b,c,d,h,w,dt
  integer::n,i,k

  !------------------------------variable globale-----------------------------------------
  w=10.0D0
  c=1.
  a=-0.6D0
  b=-a
  h=0.05

  !--------------------solution exact trace sur un intervalle ----------------------------
  call cree_tableau_pas(x,a,b,h)
  n=size(x,1)
  allocate(y(n))
  call exact_phi(w,c,x,y)   
  call ecrirefichier2d('result/exact.dat','unknown',x,y)
  
  !-----------------------------conditions initiales Euler explicite -----------------------
  M(1)=y(1)
  M(2)=exact_Dphi(w,c,x(1))   
  allocate(ye(n))
  ye(1)=M(1)

  !-------------------------------------Euler explicite -----------------------------------
  do i=2,n
     call Euler_e(2,M,M_next,h,w,c)
     ye(i)=M_next(1)
     M=M_next
  end do
  call ecrirefichier2d('result/Euler_e.dat','unknown',x,ye)


  !-----------------------------conditions initiales Euler implicite -----------------------
  M(1)=y(1)
  M(2)=exact_Dphi(w,c,x(1))   
  allocate(yei(n))
  yei(1)=M(1)

  !-------------------------------------Euler implicite -----------------------------------
  do i=2,n
     call Euler_i(2,M,M_next,h,w,c)
     yei(i)=M_next(1)
     M=M_next
  end do
  call ecrirefichier2d('result/Euler_i.dat','unknown',x,yei)
  


  !------------------------------conditions initiales RK4 --------------------------------
  M(1)=y(1)
  M(2)=exact_Dphi(w,c,x(1))   
  allocate(yk(n))
  yk(1)=M(1)

  !--------------------------------------RK4---------------------------------------------
  do i=2,n
     call RK4(2,M,M_next,h,w,c)
     yk(i)=M_next(1)
     M=M_next
  end do
  call ecrirefichier2d('result/RK4.dat','unknown',x,yk)


  !------------------------------conditions initiales Adam --------------------------------

  allocate(ya(n))

  M4(1)=y(1)
  M4(2)=exact_Dphi(w,c,x(1)) 
  M3(1)=y(2)
  M3(2)=exact_Dphi(w,c,x(2))
  M2(1)=y(3)
  M2(2)=exact_Dphi(w,c,x(3)) 
  M1(1)=y(4)
  M1(2)=exact_Dphi(w,c,x(4))
  ya(1)=M4(1)
  ya(2)=M3(1)
  ya(3)=M2(1)
  ya(4)=M1(1)

  !--------------------------------------Adam4---------------------------------------------
  do i=5,n
     call Adam4(2,M1,M2,M3,M4,M_next,h,w,c)
       ya(i)=M_next(1)
       M4=M3
       M3=M2
       M2=M1
       M1=M_next
    end do
  call ecrirefichier2d('result/Adam4.dat','unknown',x,ya)

!------------------------------conditions initiales Adam prcor--------------------------------



  M4(1)=y(1)
  M4(2)=exact_Dphi(w,c,x(1)) 
  M3(1)=y(2)
  M3(2)=exact_Dphi(w,c,x(2))
  M2(1)=y(3)
  M2(2)=exact_Dphi(w,c,x(3)) 
  M1(1)=y(4)
  M1(2)=exact_Dphi(w,c,x(4))
  ya(1)=M4(1)
  ya(2)=M3(1)
  ya(3)=M2(1)
  ya(4)=M1(1)

  !--------------------------------------Adam4 prcor---------------------------------------------
  do i=5,n
     call Adam4cor(2,M1,M2,M3,M4,M_next,h,w,c)
       ya(i)=M_next(1)
       M4=M3
       M3=M2
       M2=M1
       M1=M_next
    end do
  call ecrirefichier2d('result/Adam4cor.dat','unknown',x,ya)



!-----------------------erreur euler RK4 Adam4 norme 2 et infinie-------------------------
  
  a=-0.5
  b=-a


  call cree_tableau_pas(htab,0.000010D0,0.050D0,0.000010D0)
  n=size(htab,1)
  print*,'ok'
  allocate(e_n2(n))
  allocate(RK4_n2(n))
  allocate(e_nf(n))
  allocate(RK4_nf(n))
  allocate(adam4_n2(n))
  allocate(adam4_nf(n))
  allocate(ei_n2(n))
  allocate(ei_nf(n))
  allocate(ac_n2(n))
  allocate(ac_nf(n))

  print*,'ok' 
  do i=1,n
     call erreur_Euler_e(a,b,htab(i),w,c,e_n2(i),e_nf(i))
     call erreur_Euler_i(a,b,htab(i),w,c,ei_n2(i),ei_nf(i))
     call erreur_RK4(a,b,htab(i),w,c,RK4_n2(i),RK4_nf(i))
     call erreur_Adam4(a,b,htab(i),w,c,adam4_n2(i),adam4_nf(i))
     call erreur_Adam4cor(a,b,htab(i),w,c,ac_n2(i),ac_nf(i))
  end do

  call ecrirefichier2d('result/error_euler_norme2.dat','unknown',htab,e_n2)
  call ecrirefichier2d('result/error_RK4_norme2.dat','unknown',htab,RK4_n2)
  call ecrirefichier2d('result/error_euler_nf.dat','unknown',htab,e_nf)
  call ecrirefichier2d('result/error_RK4_nf.dat','unknown',htab,RK4_nf)
  call ecrirefichier2d('result/error_Adam4_norme2.dat','unknown',htab,adam4_n2)
  call ecrirefichier2d('result/error_Adam4_nf.dat','unknown',htab,adam4_nf)
  call ecrirefichier2d('result/error_Adam4cor_norme2.dat','unknown',htab,ac_n2)
  call ecrirefichier2d('result/error_Adam4cor_nf.dat','unknown',htab,ac_nf)
  call ecrirefichier2d('result/error_euleri_norme2.dat','unknown',htab,ei_n2)
  call ecrirefichier2d('result/error_euleri_nf.dat','unknown',htab,ei_nf)


  !------------------calcul de l'intégrale de la valeur absolue de phi^2--------------------
  
  a=-2.0D0
  b=-a
  h=0.010D0
  deallocate(x)
  call cree_tableau_pas(x,a,b,h)
  n=size(x,1)
  deallocate(y)  
  allocate(y(n))
  call exact_phi(w,c,x,y)   
  call ecrirefichier2d('result/exactint.dat','unknown',x,y)
  
  do i=1,size(x,1)
     y(i)=integrale_phi2(w,c,x(i),h)
  end do
  
  call ecrirefichier2d('result/intphi2.dat','unknown',x,y)


  !-------------------------------------calcul de u(0,x)-------------------------------------
  
  a=-2.0D0
  b=-a
  h=0.010D0
  deallocate(x)
  call cree_tableau_pas(x,a,b,h)
  n=size(x,1)
  allocate(U0(n))
  call exact_u0(w,c,x,h,U0)
  allocate(r(n))
  allocate(im(n))
  call ctoreal(U0,r,im)

  call ecrirefichier3d('result/u0.dat','unknown',r,im,x)
  
  r=abs(U0)

  call ecrirefichier2d('result/modu0.dat','unknown',x,r)

  !-------------------------------------calcul de u(t1,x)-------------------------------------
  
  dt=0.00000010D0
  
  allocate(A1(n,n))
  allocate(B1(n))
  
  call cree_matriceA(A1,U0,dt,h)
  call  cree_vecteurB (B1,U0,dt,h)
  
  allocate(L(n,n))
  allocate(U(n,n))

  call LU_decomposition(A1,L,U)

  allocate(U1(n))
  call  resolve_sysLU (L,U,B1,U1)
  
  call ctoreal(U1,r,im)
  
  call ecrirefichier3d('result/u1.dat','unknown',r,im,x)
  
  r=abs(U1)
  
  call ecrirefichier2d('result/modu1.dat','unknown',x,r)
 
  do i=1,1000
  call cree_matriceA(A1,U1,dt,h)
  call  cree_vecteurB (B1,U1,dt,h)
  call  resolve_sysLU (L,U,B1,U1)
  end do

  r=abs(U1)
  
  call ecrirefichier2d('result/modu1.dat','unknown',x,r)

  call test()

  deallocate(x)
  deallocate(y)  
  deallocate(ye)
  deallocate(yk)  
  deallocate(ya)
  deallocate(htab)
  deallocate(e_n2)
  deallocate(RK4_n2)
  deallocate(e_nf)
  deallocate(RK4_nf)
  deallocate(adam4_n2)
  deallocate(adam4_nf)
  deallocate(ei_n2)
  deallocate(ei_nf)
  deallocate(ac_n2)
  deallocate(ac_nf)



end program main
