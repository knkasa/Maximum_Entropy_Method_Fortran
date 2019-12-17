program main
implicit none

!=======================================================================
! Eliashberg code that calculates optical scattering rate in "tau.dat"
! compile w/ subroutine   sub_self2   &    sub_tau
!======================================================================

! coupling constants for Z and phi 
! spin fluctuation (sf) = mmp model,  phonon (ph) = lorentzian peak
double precision, parameter :: lamz_sf = 1.0d0  ! change these ******
double precision, parameter :: lamp_sf = 1.0d0
double precision, parameter :: lamz_ph = 0.5d0
double precision, parameter :: lamp_ph = 0.5d0

double precision, parameter :: dw = 0.005d0  ! frequency grid  ****
double precision, parameter :: maxomg = 0.5d0 ! max omega in tau 
double precision, parameter :: temp = 10.0d0  ! temperature *****
double precision, parameter :: maxomg_s = 1.0d0  ! max omega in self energy

double precision ::  pi, del, sf0, ph0
double precision, allocatable, dimension(:) :: tau, a2fz, a2fp, omg, olam
double precision, allocatable, dimension(:) :: a2fp_sf, a2fp_ph
double precision, allocatable, dimension(:) :: a2fz_sf, a2fz_ph
double complex, allocatable, dimension(:) :: zbar, phi

integer :: n, numw, numw_s, numnu, rank

!---------------------------------------------------------------------------

rank = 0
pi = 3.141592654d0

numw   = int( maxomg/dw )
numw_s = int( maxomg_s/dw )

allocate(   tau(1:numw),  a2fz(1:numw),  a2fp(1:numw),   a2fz_sf(1:numw)  )
allocate(   a2fp_sf(1:numw),  a2fz_ph(1:numw),  a2fp_ph(1:numw)   )
allocate(   zbar(-numw_s:numw_s),  phi(-numw_s:numw_s)     )
allocate(   omg(1:numw), olam(1:numw)   )

do n = 1,numw
   omg(n) = dfloat(n)*dw
end do

! -------- here you can model a2f(w) ---------------------
   sf0 = 0.041d0  ! spin fluctuation energy  
   del = 0.005d0  ! broadening of phonon peak
   ph0 = 0.043d0  ! position of phonon peak

   a2fz_sf = omg/sf0/( omg*omg + sf0*sf0 )  ! mmp model for sf
   a2fp_sf = omg/sf0/( omg*omg + sf0*sf0 )

   a2fz_ph = del/pi/( del*del + (omg-ph0)**2 )  ! lorentzian for ph
   a2fp_ph = del/pi/( del*del + (omg-ph0)**2 ) 
!----------------------------------------------------------

   ! normalizing
   a2fz_sf = lamz_sf * a2fz_sf/2.0d0/dw/sum( a2fz_sf/omg )
   a2fp_sf = lamp_sf * a2fp_sf/2.0d0/dw/sum( a2fp_sf/omg )
   a2fz_ph = lamz_ph * a2fz_ph/2.0d0/dw/sum( a2fz_ph/omg )
   a2fp_ph = lamp_ph * a2fp_ph/2.0d0/dw/sum( a2fp_ph/omg )

   a2fz = a2fz_sf + a2fz_ph
   a2fp = a2fp_sf + a2fp_ph

!open(11,file='a2f.dat')
!open(11,file='a2f-check.dat')
!do n = 1,numw
!    read(11,*)   omg(n), a2fz(n), a2fp(n)
!   write(11,*)  omg(n), a2fz(n), a2fp(n)
!end do
!close(11)
  
  call sub_self2( rank, a2fz, a2fp, numw, numw_s, maxomg_s, temp, zbar, phi )
  call sub_tau2(   zbar, phi, numw_s, maxomg_s, numw, temp, tau, olam  ) 

 if( maxval( real(phi) )<0.001d0 ) print*,' Note: system not superconducting'
 

open(21,file='tau-check.dat')
do n = 1,numw
  write(21,*)  omg(n), tau(n), olam(n)  
end do
close(21)

print*, '  ****** Done ******   '

end program main


