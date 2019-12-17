program cond1
implicit none


!========================================================================
! Tc calculation for infinite band model. d-wave assumed. 
!========================================================================

!-----------------------------------------------------------------------
!   Run 'a2f.f95' to obtain spectral density 
!-----------------------------------------------------------------------
 
integer, parameter :: numv = 200  ! # of freq in a2f (see mycond "numnu") *****
double precision, parameter :: vmax = 0.4d0  !*********
 ! max freq in spectral density (see nuscale)

integer, parameter :: numt = 80  ! # of temp points  ****
double precision, parameter :: ti = 20.0d0  ! initial temp for plot *****
double precision, parameter :: tf = 200.0d0 ! final temp for plot   *****


integer, parameter :: ntheta = 32  ! # of theta points in theta integral
double precision, parameter :: ztophi=2.0 ! thiss should be 2.0
integer :: maxit = 500  ! # of max iteration

double precision ::  dv, dw,  beta, pi, dtheta, phi, exc, diff
double precision ::  bs, fm, fp, nu, teta, w, kb
integer :: i,j,l,m ,t, nmax, iter
double precision, dimension(1:numv) :: v, Bz, Bph
!double precision, dimension(-nmax:nmax-1) :: wn
!double precision, dimension(-nmax:nmax) :: wm, ilamz, ilamph 
double precision, dimension(0:ntheta) :: theta, fac
double complex, parameter :: ii = cmplx(0.0,1.0d0)
double complex :: exc2, lamzm, lamzp, lampm, lampp
!double complex, dimension(0:nmax-1) :: iwm
double complex, allocatable, dimension(:) :: z, gap, zp, gapp, iwm
double precision, allocatable, dimension(:) :: wm, ilamz, ilamph, wn

double precision, parameter :: sig = 1.0e-7 ! convergence criterion
double precision  ::  temp, delta

!-------------------------------------------------------------------------

pi = 3.141592654d0
kb = 8.617e-5
dv = vmax/dfloat(numv)

open(21,file='a2f.OUT')  
!open(21,file='a2f2.dat') ! input file for a2f(w) *********
do m=1,numv
        read(21,*) v(m), Bz(m), Bph(m)  !spectral density
end do
close(21)

do m = 0,ntheta
        dtheta=0.25d0*pi/dfloat(ntheta)
        theta(m) = dtheta*dfloat(m)   ! theta
end do

do m = 0,ntheta,2  ! the factor that appears in simpson's rule
        fac(m) = 4.0d0/3.0d0
        fac(m+1) = 2.0d0/3.0d0
end do
fac(0) = 1.0d0/3.0d0
fac(ntheta) = 1.0d0/3.0d0



!----------  Initializing  -------------------------------------------

open(13,file='gap.dat')
do t = 0,numt
    temp = dfloat(t)*(tf-ti)/dfloat(numt) + ti
    beta = 1.0d0/kb/temp
    nmax = int( vmax*beta/pi ) 
    print*, ' '
    print*, ' # of matsubara freq = ', nmax

allocate( iwm(0:nmax-1), wn(-nmax:nmax-1), wm(-nmax:nmax)          )
allocate( ilamz(-nmax:nmax), ilamph(-nmax:nmax), z(-nmax:nmax-1)   )
allocate( gap(-nmax:nmax-1), zp(-nmax:nmax-1), gapp(-nmax:nmax-1)  )


! ( preparing wn, Bsf, omg, iwm )
do m = -nmax,nmax-1
        wn(m) = pi*(2.0d0*m+1.0d0)/beta ! wn = fermion matsubara poles
end do
do m = 0,nmax-1  !  iwm = imaginary  matsubara sum in 1st term
        iwm(m) = dcmplx( 0.0d0, wn(m)  ) 
end do
do m = -nmax,nmax
        wm(m) = pi*2.0d0*m/beta  ! boson matsubara poles
end do


! calculating lambda(i*wm) 
ilamz = 0.0d0
ilamph= 0.0d0
do i = -nmax,nmax  ! need boson matsubara poles
do j = 1,numv
    ilamz(i)  = ilamz(i) + 2.0d0*dv*v(j)*Bz(j)/( v(j)*v(j) + wm(i)*wm(i) )
    ilamph(i) = ilamph(i)+ 2.0d0*dv*v(j)*Bph(j)/( v(j)*v(j) + wm(i)*wm(i) )
end do       
end do


!-------------  first getting imaginary self energy  --------------------
! ( obtaining self energy with imaginary matsubara freq dependence )


do m = -nmax,nmax-1  
        gap(m) = dcmplx(0.1d0,0.0) ! initial guess of gap & Z
          z(m) = dcmplx(wn(m),0.0)
end do

iter = 1
diff = 1.0d0
do  while  ( diff > sig .and. iter < maxit )     ! starting iteration
     
       gapp = gap
         zp = z
do i = -nmax,nmax-1
          z(i) = wn(i)
        gap(i) = dcmplx( 0.0d0, 0.0d0 )
do j = -nmax,nmax-1
       if ( abs(j-i) <= nmax ) then
    do l = 0,ntheta

        phi = cos(2.0d0*theta(l))
        exc = real( sqrt( zp(j)*zp(j) + gapp(j)*gapp(j)*phi*phi ) )
   
      z(i) = z(i) +  pi/beta*ilamz(j-i)*zp(j)*dtheta/(exc*0.25d0*pi)*fac(l) 
      gap(i) = gap(i) + pi*ztophi/beta*ilamph(j-i)*gapp(j)*phi**2*dtheta &
           /(exc*0.25d0*pi)  * fac(l)

    end do    ! end l
       end if
end do  ! end j
end do  ! end i
        diff = maxval( abs( real( z - zp )  )  )
        print*, iter,  diff
        iter = iter+1 
end do  ! end while
     
    
     delta = real( gap(0) )/real( Z(0) ) *pi/beta
     write(13,*) temp, delta
     print*, '------------------------------------------------- '
     print*, ' temp and gap = ', temp, delta
     print*, '------------------------------------------------ '

deallocate( iwm, wn, wm   )
deallocate( ilamz, ilamph, z   )
deallocate( gap, zp, gapp  )


end do  ! end temp loop
close(13)
!-----------------  end of imaginary self energy  -------------------


end  program cond1 

