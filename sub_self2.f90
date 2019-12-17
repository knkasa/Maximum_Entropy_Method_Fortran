subroutine  sub_self2( rank, a2fz, a2fp, numv, numw, wmax, temp,  z, gap )
implicit none


!========================================================================
! This code calculates self energy with real freq dependence
! d-wave single band cuparates is assumed  
! See Equation(16), Carbotte PRB 71, 054506 (2005)
! For optical conductivities, see cond1.f95 
!=======================================================================


integer :: maxit = 200  ! # of max iteration
double precision, parameter :: conv = 1.0e-7  ! convergence criteria

integer, parameter :: nmax= 100    ! max N in matsubara frequency 
integer, parameter :: ntheta = 16  ! # of theta points in theta integral

double precision ::  dv, dw, T, beta, pi, dtheta, phi, exc, diff
double precision ::  bs, fm, fp, nu, teta, w, ztophi, wmax, temp
integer :: i,j,l,m, iter, numw, numv, rank

double precision, dimension(1:numv) :: v, Bz, Bph, a2fz, a2fp
double precision, dimension(-nmax:nmax-1) :: wn
double precision, dimension(-nmax:nmax) :: wm, ilamz, ilamph 
double precision, dimension(0:ntheta) :: theta, fac
double precision, dimension(-numw:numw) :: omg,  dos
double complex, parameter :: ii = cmplx(0.0,1.0d0)
double complex :: exc2, lamzm, lamzp, lampm, lampp
double complex, dimension(0:nmax-1) :: iwm
double complex, dimension(-numw:numw) :: cz, cp

double complex, dimension(-nmax:nmax-1) :: iz, igap, iz2, igap2
double complex, dimension(-numw:numw) :: z, gap, zp, gapp


T = 8.617E-5*temp  ! Temperature T=10 degree ******** 
beta = 1.0d0/T
pi = 3.141592654d0
ztophi = 2.0d0   !  z to phi ratio (it needs to be two)
dtheta=0.25d0*pi/dfloat(ntheta)


!--------------  preparing wn, theta, Bsf, omg, iwm ---------------

do m = -nmax,nmax-1
        wn(m) = pi*(2.0d0*m+1.0d0)/beta ! wn = fermion matsubara poles
end do

do m = 0,ntheta
        theta(m) = dtheta*dfloat(m)   ! theta
end do


Bz = a2fz
Bph = a2fp         
dw = wmax/dfloat(numw) 
dv = dw


do m = -numw,numw
        omg(m) = dw*dfloat(m)   ! omg =  real freuency for final result
end do

do m = 0,nmax-1  !  iwm = imaginary  matsubara sum in 1st term of zbar
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
    nu = dv*dfloat(j)
    ilamz(i)  = ilamz(i) + 2.0d0*dv*nu*Bz(j)/( nu*nu + wm(i)*wm(i) )
    ilamph(i) = ilamph(i)+ 2.0d0*dv*nu*Bph(j)/( nu*nu + wm(i)*wm(i) )
end do       
end do


!-------------  first getting imaginary self energy  --------------------
! ( obtaining self energy with imaginary matsubara freq dependence )

!print*, ' '
!print*, 'Obtaining self energy on imaginary axis'

!if(rank==0)then
!print*, ' '
!print*, 'Calculatin Imaginary part'
!print*, ' '
!end if
do m = -nmax,nmax-1  
        igap(m) = dcmplx(0.01d0,0.0) ! initial guess of gap & Z
        iz(m) = dcmplx(wn(m),0.0)
end do



! starting iteration
diff = 1.0d0
iter = 1
do while( diff>conv .and. iter<maxit )
        igap2 = igap
        iz2 = iz
do i = -nmax,nmax-1
        iz(i) = wn(i)
        igap(i) = dcmplx( 0.0d0, 0.0d0 )
do j = -nmax,nmax-1
       if ( abs(j-i) <= nmax ) then
    do l = 0,ntheta

        phi = cos(2.0d0*theta(l))
        exc = real( sqrt( iz2(j)*iz2(j) + igap2(j)*igap2(j)*phi*phi ) )
   
     iz(i) = iz(i) +  pi*T*ilamz(j-i)*iz2(j)*dtheta/(exc*0.25d0*pi) *fac(l) 
      igap(i) = igap(i) + pi*ztophi*T*ilamph(j-i)*igap2(j)*phi**2*dtheta &
           /(exc*0.25d0*pi) * fac(l)

     end do    ! end l
       end if
end do  ! end j
end do  ! end i
        diff = abs( maxval( real( iz-iz2  ) )  )
    !if(rank==0)then    
    ! print*,   iter, diff
    !endif
        iter  = iter+1

   if(rank==0)then
     if(iter==200) print*, 'Imaginary part not converging', iter, diff
   end if

end do  ! end while


!open(11,file='igap2.OUT')
!open(12,file='iz2.OUT')
!do m = -nmax,nmax-1
!        write(11,*) wn(m), dreal(igap(m)), dimag(igap(m))
!        write(12,*) wn(m), dreal(iz(m)), dimag(iz(m))
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
   nu = dfloat(l)*dv
  lamzm=lamzm+Bz(l)*dv*(1.0d0/(omg(i)-iwm(j)+nu) -1.0d0/(omg(i)-iwm(j)-nu))
  lamzp=lamzp+Bz(l)*dv*(1.0d0/(omg(i)+iwm(j)+nu) -1.0d0/(omg(i)+iwm(j)-nu))
 lampm=lampm+Bph(l)*dv*(1.0d0/(omg(i)-iwm(j)+nu) -1.0d0/(omg(i)-iwm(j)-nu))
 lampp=lampp+Bph(l)*dv*(1.0d0/(omg(i)+iwm(j)+nu) -1.0d0/(omg(i)+iwm(j)-nu))

end do  
do l = 0,ntheta
        phi = cos(2.0d0*theta(l))
        exc = dreal( sqrt( iz(j)*iz(j) + igap(j)*igap(j)*phi*phi ) ) 
     
  cz(i) = cz(i) +  ii*pi*T*(lamzm-lamzp)*iz(j)*dtheta / (exc*0.25d0*pi)*fac(l)
  cp(i)= cp(i)+ ztophi*pi*T*(lampm+lampp)*igap(j)*phi**2*dtheta & 
           /(exc*0.25d0*pi) * fac(l)

end do    ! end l
end do    ! end j
end do    ! end i

!open(31,file='cp2.OUT')
!open(32,file='cz2.OUT')
!do m = -numw,numw
!        write(31,*)  omg(m),  dreal(cp(m)), dimag(cp(m))
!        write(32,*) omg(m), dreal(cz(m)), dimag(cz(m))
!end do  
!close(31)
!close(32)

!-----------------  end of calculating 1st term  -----------------------


!---------------  calculating 2nd term  -----------------------------------

!print*, ' '
!print*, 'Obtaining self energy on real axis'

   z = cz+0.01d0*ii    ! addin small imaginary 
 gap = cp+0.01d0*ii  ! (this will make convergence faster)

! starting iteration
diff=1.0d0
iter=1

!if(rank==0)then
!print*, ' '
!print*, 'Calculating real part'
!print*, ' '
!end if
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
   z(i) = z(i) +  ii*pi*dtheta*(bs+fp)*Bz(j)*dv*zp(i+j)/(0.25d0*pi*exc2)*fac(l)
     gap(i) = gap(i) + ztophi*ii*dtheta*pi*(bs+fp)*Bph(j)*dv*gapp(i+j)* &
           phi*phi /(0.25d0*pi*exc2) * fac(l)
          end if

          if( abs(i-j) <= numw ) then
       exc2 = sqrt( zp(i-j)*zp(i-j) - gapp(i-j)*gapp(i-j)*phi*phi )
          if( dimag(exc2) < 0.0d0 ) exc2 = -exc2

   z(i)  = z(i) + ii*pi*dtheta*(bs+fm)*Bz(j)*dv*zp(i-j)/(0.25d0*pi*exc2)*fac(l)
   gap(i) = gap(i)+ ztophi*ii*dtheta*pi*(bs+fm)*Bph(j)*dv*gapp(i-j)* &
           phi*phi /(0.25d0*pi*exc2) * fac(l)
          end if

end do   ! end l
end do   ! end j
end do   ! end i
        diff = abs( maxval( dimag(z-zp) )  )
     ! if(rank==0)then
     !  print*, iter,  diff
     ! endif
        iter = iter+1

   if(rank==0)then
     if(iter==200) print*, 'Real part not converging', iter, diff
   end if

end do   ! end while

!open(41,file='phi2.OUT')
!open(42,file='omgbar2.OUT')
!do m = -numw,numw
!        write(41,*) omg(m), dreal(gap(m)), dimag(gap(m))
!        write(42,*) omg(m), dreal(z(m)), dimag(z(m))
!end do
!close(41)
!close(42)

!------------  end of 2nd term  --------------------------------------



!------------- Producing output files for gap(w) & Z(w) -----------------

!open(51,file='realz.OUT')
!open(52,file='rgap.OUT')
!do m = -numw,numw
!        w = dfloat(m)*dw
!      if ( m /= 0 ) then
!        z(m) = z(m)/w
!        gap(m) = gap(m)/z(m)
!      end if
!      if ( m == 0 ) then
!        z(m) = 0.0d0
!        gap(m) = gap(1)
!      end if
!write(51,*) w, dreal(z(m)), dimag(z(m))
!write(52,*)  w, w, real(gap(m)) 
!end do
!close(51)
!close(52)


!------------ calculating density of states ------------------------------

!open(61,file='dos.OUT')
!do m = 0,numw
!        dos(m) = 0.0d0
!        w = dfloat(m)*dw
!        dtheta =  0.5d0*pi/dfloat(ntheta)
!do l = 0,ntheta
!        teta = dfloat(l)*dtheta
!        phi = cos(2.0d0*teta)
        
!        exc2 = sqrt( w*w - gap(m)*gap(m)*phi*phi )
!        dos(m) = dos(m) + dreal(w*dtheta/exc2)/(0.5d0*pi)
!end do  ! end l
!        write(61,*) w, dos(m)
!end do  ! end m
!close(61)

return
end  subroutine  sub_self2 





