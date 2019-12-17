program cond1
implicit none


!========================================================================
! This code calculates self energy with real freq dependence
! d-wave single band cuparates is assumed  
! See Equation(16), Carbotte PRB 71, 054506 (2005)
! For optical conductivities, see cond1.f95 
!=======================================================================

!-----------------------------------------------------------------------
!   Run 'a2f.f95' to obtain spectral density 
!----------------------------------------------------------------------
 
integer :: maxit = 500  ! # of max iteration
double precision, parameter :: conv = 1.0e-7  ! convergence criteria

integer, parameter :: numv = 400  ! # of v in a2f(v)(see a2f.f90 "numnu")***
double precision, parameter :: vmax = 0.4d0  ! see a2f.f90 ******
 ! ****max freq in spectral density (see nuscale)

integer, parameter :: nmax = 100    ! max N in matsubara frequency
integer, parameter :: numw = 1000   !  # of w points in real axis (10*numv)*** 
double precision, parameter :: wmax = 1.0d0  ! 4.0  max omega in real axis
integer, parameter :: ntheta = 32  ! # of theta points in theta integral


double precision ::  dv, dw, T, beta, pi, dtheta, phi, exc, diff
double precision ::  bs, fm, fp, nu, teta, w, ztophi
integer :: i,j,l,m, iter
double precision, dimension(1:numv) :: v, Bz, Bph
double precision, dimension(-nmax:nmax-1) :: wn
double precision, dimension(-nmax:nmax) :: wm, ilamz, ilamph 
double precision, dimension(0:ntheta) :: theta, fac
double precision, dimension(-numw:numw) :: omg,  dos
double complex, parameter :: ii = cmplx(0.0,1.0d0)
double complex :: exc2, lamzm, lamzp, lampm, lampp
double complex, dimension(0:nmax-1) :: iwm
double complex, dimension(-numw:numw) :: cz, cp
double complex, allocatable, dimension(:) :: z, gap, zp, gapp

T = 8.617E-5*10.0d0  ! Temperature T=10 degree ******** 
beta = 1.0d0/T
pi = 3.141592654d0
ztophi = 2.0d0   !  z to phi ratio (it needs to be two)

!--------------  preparing wn, theta, Bsf, omg, iwm ---------------

do m = -nmax,nmax-1
        wn(m) = pi*(2.0d0*m+1.0d0)/beta ! wn = fermion matsubara poles
end do

do m = 0,ntheta
        dtheta=0.25d0*pi/dfloat(ntheta)
        theta(m) = dtheta*dfloat(m)   ! theta
end do

! open(21,file='a2f-model.dat') ! input data a2f(w) **********
  open(21,file='a2f.OUT')
do m=1,numv
        read(21,*) v(m), Bz(m), Bph(m)  !spectral density
end do
close(21)

dv = vmax/dfloat(numv) 
dw = wmax/dfloat(numw)  !  *** dv=dw *** should be true
print*, ' '
print*, 'make sure dv=dw(see code)'
print*, 'dv=', dv, 'dw=', dw
dv = dw

do m = -numw,numw
        omg(m) = dw*dfloat(m)   ! omg =  real freuency for final result
end do

do m = 0,nmax-1  !  iwm = imaginary  matsubara sum in 1st term
        iwm(m) = dcmplx( 0.0d0, wn(m)  ) 
end do

do m = -nmax,nmax
        wm(m) = pi*2.0d0*m/beta  ! boson matsubara poles
end do

do m = 0,ntheta,2  ! the factor that appears in simpson's rule
        fac(m) = 4.0d0/3.0d0
        fac(m+1) = 2.0d0/3.0d0
end do
fac(0) = 1.0d0/3.0d0
fac(ntheta) = 1.0d0/3.0d0

   !fac = 1.0d0

!----------  calculating lambda(i*wm)  ------------------------------

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

print*, ' '
print*, 'Obtaining self energy on imaginary axis'

allocate( gap(-nmax:nmax-1), z(-nmax:nmax-1) )
allocate( gapp(-nmax:nmax-1), zp(-nmax:nmax-1) )

do m = -nmax,nmax-1  
        gap(m) = dcmplx(0.01d0,0.0) ! initial guess of gap & Z
        z(m) = dcmplx(wn(m),0.0)
end do


! starting iteration
diff = 1.0d0
iter = 1
do while( diff>conv .and. iter<maxit )
        gapp = gap
        zp = z
do i = -nmax,nmax-1
        z(i) = wn(i)
        gap(i) = dcmplx( 0.0d0, 0.0d0 )
do j = -nmax,nmax-1
       if ( abs(j-i) <= nmax ) then
       do l = 0,ntheta

        phi = cos(2.0d0*theta(l))
        exc = dreal( sqrt( zp(j)*zp(j) + gapp(j)*gapp(j)*phi*phi ) )
   
     z(i) = z(i) +  pi*T*ilamz(j-i)*zp(j)*dtheta/(exc*0.25d0*pi) *fac(l)
      gap(i) = gap(i) + pi*ztophi*T*ilamph(j-i)*gapp(j)*phi**2*dtheta &
           /(exc*0.25d0*pi) * fac(l)

       end do    ! end l
       end if
end do  ! end j
end do  ! end i
        diff = abs( dreal(gap(0) - gapp(0))  )
        print*,   iter, diff
        iter  = iter+1
end do  ! end while


!open(11,file='igap.OUT')
!open(12,file='iz.OUT')
!do m = -nmax,nmax-1
!        write(11,*) wn(m), dreal(gap(m)), dimag(gap(m))
!        write(12,*) wn(m), dreal(z(m)), dimag(z(m))
!end do 
!close(11)
!close(12)

!-----------------  end of imaginary self energy  -------------------




!------------------  calculating 1st term  ---------------------------

do i = -numw,numw  ! real frequency omg
        cz(i) = dcmplx( omg(i), 0.0d0 )
        cp(i) = ( 0.0d0, 0.0d0)
do j = 0,nmax-1    ! summing over matsubara poles
        lamzm = dcmplx( 0.0d0, 0.0d0) ! lamda z minus
        lamzp = dcmplx( 0.0d0, 0.0d0) ! lamda z plus
        lampm = dcmplx( 0.0d0, 0.0d0) ! lamda phi minus
        lampp = dcmplx( 0.0d0, 0.0d0) ! lamda phi plus
do l = 1,numv

  lamzm=lamzm+Bz(l)*dv*(1.0d0/(omg(i)-iwm(j)+v(l)) -1.0d0/(omg(i)-iwm(j)-v(l)))
  lamzp=lamzp+Bz(l)*dv*(1.0d0/(omg(i)+iwm(j)+v(l)) -1.0d0/(omg(i)+iwm(j)-v(l)))
 lampm=lampm+Bph(l)*dv*(1.0d0/(omg(i)-iwm(j)+v(l)) -1.0d0/(omg(i)-iwm(j)-v(l)))
 lampp=lampp+Bph(l)*dv*(1.0d0/(omg(i)+iwm(j)+v(l)) -1.0d0/(omg(i)+iwm(j)-v(l)))

end do  
do l = 0,ntheta
        phi = cos(2.0d0*theta(l))
        exc = dreal( sqrt( z(j)*z(j) + gap(j)*gap(j)*phi*phi ) ) 
     
    cz(i) = cz(i) + ii*pi*T*(lamzm-lamzp)*z(j)*dtheta/(exc*0.25d0*pi) *fac(l)
  cp(i) = cp(i)+ ztophi*pi*T*(lampm+lampp)*gap(j)*phi**2*dtheta &
          /(exc*0.25d0*pi) *fac(l)

end do    ! end l
end do    ! end j
end do    ! end i

!open(31,file='cp.OUT')
!open(32,file='cz.OUT')
!do m = -numw,numw
!        write(31,*)  omg(m),  dreal(cp(m)), dimag(cp(m))
!        write(32,*) omg(m), dreal(cz(m)), dimag(cz(m))
!end do  
!close(31)
!close(32)

!-----------------  end of calculating 1st term  -----------------------




!---------------  calculating 2nd term  -----------------------------------

print*, ' '
print*, 'Obtaining self energy on real axis'

deallocate( z, gap, zp, gapp )
allocate( z(-numw:numw), gap(-numw:numw), zp(-numw:numw), gapp(-numw:numw) )

 z = cz+0.01d0*ii    ! addin small imaginary 
 gap = cp+0.01d0*ii  ! (this will make convergence faster)

! starting iteration
diff=1.0d0
iter=1
do while( diff>conv .and. iter<maxit )
        zp = z
        gapp = gap
do i = -numw,numw
        w = dfloat(i)*dw        
        z(i) = cz(i)
        gap(i) = cp(i)
do j = 1,numv
        nu = dfloat(j)*dw 

        bs = 1.0d0/( dexp(beta*nu) - 1.0d0 )  ! bose function
        fm = 1.0d0/( dexp(beta*(nu-w)) + 1.d00 )  ! fermi func f(v-wm)
        fp = 1.0d0/( dexp(beta*(nu+w)) + 1.0d0 )  ! fermi func f(v+wm)

do l = 0,ntheta
        teta = dfloat(l)*dtheta
        phi = cos(2.d00*teta)
        
          if ( abs(i+j) <= numw ) then
       exc2 = sqrt( zp(i+j)*zp(i+j) - gapp(i+j)*gapp(i+j)*phi*phi )
           if ( dimag(exc2) < 0.0d0 ) exc2 = -exc2 
            ! (above) this way, no need for end if
   z(i) = z(i) + ii*pi*dtheta*(bs+fp)*Bz(j)*dv*zp(i+j)/(0.25d0*pi*exc2)*fac(l)
     gap(i) = gap(i) + ztophi*ii*dtheta*pi*(bs+fp)*Bph(j)*dv*gapp(i+j)* &
           phi*phi /(0.25d0*pi*exc2) *fac(l)
          end if

          if( abs(i-j) <= numw ) then
       exc2 = sqrt( zp(i-j)*zp(i-j) - gapp(i-j)*gapp(i-j)*phi*phi )
          if( dimag(exc2) < 0.0d0 ) exc2 = -exc2

   z(i)  = z(i) + ii*pi*dtheta*(bs+fm)*Bz(j)*dv*zp(i-j)/(0.25d0*pi*exc2)*fac(l)
    gap(i) = gap(i)+ ztophi*ii*dtheta*pi*(bs+fm)*Bph(j)*dv*gapp(i-j)* &
           phi*phi /(0.25d0*pi*exc2) *fac(l)
          end if

end do   ! end l
end do   ! end j
end do   ! end i
        diff = abs( dimag(gap(-100)-gapp(-100) ) )
        print*, iter,  diff
        iter = iter+1
end do   ! end while

open(41,file='phi.OUT')
open(42,file='omgbar.OUT')
do m = -numw,numw
        write(41,*) omg(m), dreal(gap(m)), dimag(gap(m))
        write(42,*) omg(m), dreal(z(m)), dimag(z(m))
end do
close(41)
close(42)

!------------  end of 2nd term  --------------------------------------



!------------- Producing output files for gap(w) & Z(w) -----------------

!open(51,file='realz.OUT')
open(52,file='rgap.OUT')
do m = -numw,numw
        w = dfloat(m)*dw
      if ( m /= 0 ) then
        z(m) = z(m)/w
        gap(m) = gap(m)/z(m)
      end if
      if ( m == 0 ) then
        z(m) = 0.0d0
        gap(m) = gap(1)
      end if
!write(51,*) w, dreal(z(m)), dimag(z(m))
write(52,*)  w, w, real(gap(m)) 
end do
!close(51)
close(52)


!------------ calculating density of states ------------------------------

open(61,file='dos.OUT')
do m = 0,numw
        dos(m) = 0.0d0
        w = dfloat(m)*dw
        dtheta =  0.5d0*pi/dfloat(ntheta)
do l = 0,ntheta
        teta = dfloat(l)*dtheta
        phi = cos(2.0d0*teta)
        
        exc2 = sqrt( w*w - gap(m)*gap(m)*phi*phi )
        dos(m) = dos(m) + dreal(w*dtheta/exc2)/(0.5d0*pi)
end do  ! end l
        write(61,*) w, dos(m)
end do  ! end m
close(61)


end  program cond1 





