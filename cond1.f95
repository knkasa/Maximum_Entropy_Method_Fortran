program cond
implicit none

!-----------------------------------------------------------------------------
!  This will calculate optical conductivity 
!  See eq(15) Carbotte, PRB 71 054506
!  Self energies are calculated from model05.f95
!----------------------------------------------------------------------------

integer, parameter :: numw = 1000   ! max # of omega (see model05.f95) ***
integer, parameter :: ntheta = 32  ! # of theta points
double precision, parameter :: wmax = 1.0d0 ! omega max (see model05.f95)***

double precision :: pi, dtheta, theta, phi, T, beta, dw
double precision :: lamop, tauop
double precision :: mag, omg, nu, z1, z2, gap1, gap2
double precision, dimension(0:ntheta) :: fac

double complex :: Jp, Jm, den 
double complex :: ctmp, frac, cone, ii
double complex :: E0, Ep, Em, NN, Np, Nm, P, Pp, Pm
double complex, dimension(-numw:numw) :: z, gap
double complex, dimension(0:numw) :: sgm
integer :: i, j, l, m


T = 10.0d0*0.00008617d0  ! k*T  T=10 degree (see mycond.f95)************
pi = 3.141592654d0
beta = 1.0d0/T
dw = wmax/dfloat(numw)  
cone = dcmplx(1.0d0, 0.0d0)
ii = dcmplx(0.0d0, 1.0d0)
dtheta = 0.25d0*pi/dfloat(ntheta)


do m = 0,ntheta,2  ! the factor that appears in simpson's rule
        fac(m) = 4.0d0*dtheta/3.0d0
        fac(m+1) = 2.0d0*dtheta/3.0d0
end do

fac(0) = dtheta/3.0d0
fac(ntheta) = dtheta/3.0d0

!fac = dtheta

open(11,file='omgbar.OUT')
open(12,file='phi.OUT')
do m = -numw,numw
        read(11,*) omg, z1, z2
        read(12,*) omg, gap1, gap2
        z(m) = dcmplx( z1, z2 )
        gap(m) = dcmplx( gap1, gap2 )
end do
close(11)
close(12)



do i = 1,numw/2   !numw/4  ! omega for plotting**************
        sgm(i) = dcmplx(0.0d0, 0.0d0)
        nu = dfloat(i)*dw
do j = 0,numw    ! omega for integration
        omg = dfloat(j)*dw
        ctmp = dcmplx(0.0d0, 0.0d0)
do l = 0,ntheta  ! theta for integration
        theta = dfloat(l)*dtheta
        phi = cos(2.0d0*theta)


        E0 = sqrt( z(j)*z(j) - gap(j)*gap(j)*phi*phi )
      if (dimag(E0) < 0.0d0 )   E0 = -E0        

         NN = z(j)/E0
         P = gap(j)*phi/E0

        
    if ( abs(j+i) <= numw ) then
        Ep = sqrt( z(j+i)*z(j+i) - gap(j+i)*gap(j+i)*phi*phi )
       if ( dimag(Ep) < 0.0d0 ) Ep = -Ep 
        Np = z(j+i)/Ep
        Pp = gap(j+i)*phi/Ep
        
        den = E0+Ep
        mag = dreal(den*dconjg(den))
       if ( mag /= 0.0d0) Jp = (cone - NN*Np - P*Pp)/den
        den = dconjg(E0) - Ep
        mag = dreal(den*dconjg(den))
       if (mag /= 0.0d0) Jp=Jp+(cone+dconjg(NN)*Np+dconjg(P)*Pp)/den
    else
        Jp = 0.0d0
    end if


    if( abs(j-i) <= numw )  then
        Em = sqrt( z(j-i)*z(j-i) - gap(j-i)*gap(j-i)*phi*phi )
       if ( dimag(Em) < 0.0d0 ) Em = -Em
        Nm = Z(j-i)/Em
        Pm = gap(j-i)*phi/Em    

        den = E0+Em
        mag = dreal(den*dconjg(den))
       if (mag /= 0.0d0) Jm = (cone - NN*Nm - P*Pm)/den
        den = dconjg(E0)-Em
        mag = dreal(den*dconjg(den))
       if (mag /= 0.0d0) Jm=Jm+(cone+dconjg(NN)*Nm+dconjg(P)*Pm)/den
    else
        Jm = 0.0d0
    end if

        ctmp = ctmp + 0.5d0*(Jp+dconjg(Jm))*fac(l)    

end do   ! end l (theta integral)
        
        sgm(i) = sgm(i) + dtanh(beta*omg*0.5d0)*ctmp*dw/(0.25d0*pi)

end do   ! end j (omega integral)
        sgm(i) = sgm(i)*ii/nu
end do   ! end i 


open(21,file='sigma.OUT')
do m = 1,numw/2    !numw/4   ! this needs to be consistent with above  *******
        nu = dfloat(m)*dw
      
        lamop = dimag(cone/sgm(m)) + nu
        tauop = dreal(cone/sgm(m))
        
write(21,*) nu, dreal(sgm(m)), dimag(sgm(m)), lamop, tauop
end do
close(21)



end program cond
