program main
implicit none
include 'mpif.h'

!===================================================================
! Calculating kernel using full Eliashberg equation.
! Temperature and coupling constant(lambda) are input parameters.
! This assumes kernel is independent of a2f(w)
! Compile with subroutine mpif90 kernel.f90 sub_self.f90 sub_tau.f90
! Midpoint averaging is used to calculate derivative.
!==================================================================== 

! Note:  Kernel is dependent on temperature and coupling strength(lambda)
! kernel = K(T,lambda)  Set lambda appropriately to match Tc with experiment

integer :: rank, numpc, ierr  ! these for mpi

double precision, parameter ::  dw = 0.001d0  ! freqency grid  *******
double precision, parameter :: maxomg = 0.5d0 ! max energy for kernel *****
double precision, parameter :: temp = 10.0d0  ! temperature ********
double precision, parameter :: lambda = 1.54d0 ! coupling constant
! 1.252d0

double precision, parameter :: diff = 0.0001d0 ! small deviation in derivative
double precision, parameter :: maxomg_s = 1.0d0  ! max energy in self energy   

double precision :: ph0, ph_with, pi, kb, beta, nu
double precision, allocatable, dimension(:) :: omg, rho, rho_sf, rho_ph
double precision, allocatable, dimension(:) :: tau, tau2, dtau, tau1
double complex, allocatable, dimension(:) :: zbar,  phi
double precision, allocatable, dimension(:,:) :: rho_mat1, kernel, rho_mat2

integer :: n, i, j, numw, numw_s


!-------------------------------------------------------------------

call mpi_init(ierr)
call mpi_comm_rank(mpi_comm_world, rank, ierr)
call mpi_comm_size(mpi_comm_world, numpc, ierr)

kb = 8.617e-5
beta = 1.0d0/kb/temp
pi =  3.141592654d0
ph0 = 0.041d0  ! position of phonon peak in a2f(w)
ph_with = 0.005d0  ! width of phonon peak in a2f(w)

numw = int( maxomg/dw )
numw_s = int( maxomg_s/dw )

allocate( omg(1:numw), rho(1:numw), rho_sf(1:numw), rho_ph(1:numw)  )
allocate( tau(1:numw), tau2(1:numw), dtau(1:numw), tau1(1:numw)     )
allocate( zbar(-numw_s:numw_s),  phi(-numw_s:numw_s)      )
allocate( rho_mat1(1:numw,1:numw), kernel(1:numw,1:numw)   )
allocate( rho_mat2(1:numw,1:numw)  )


do n = 1,numw
   omg(n)  = dfloat(n)*dw
end do


! initializing a2f(w) = rho(w)
rho_sf = omg/ph0/( omg*omg + ph0*ph0 )
rho_ph = ph_with/pi/( ph_with**2 + (omg-ph0)**2 )


!do n = 1,numw
 !  nu = dfloat(n)*dw
 !  if ( nu>0.4000001d0 ) then  ! zero-ing a2f(w) just for model
    ! rho_sf(n) = 0.0d0
    ! rho_ph(n) = 0.0d0
!   end if
!end do
 
! normalizing rho(w) (phonon contribution should be zero for better result)
rho = lambda * rho_sf/2.0d0/dw/sum(rho_sf/omg)  !  & 
                     !      +  0.0d0 * rho_ph/2.0d0/dw/sum(rho_ph/omg)

deallocate( rho_sf, rho_ph, omg )

!open(21,file='test.dat')
!do n = 1,numw
!   write(21,*) n, rho(n)
!end do
!close(21)


! Note the order is important!!  
!call sub_self( rho, numw, numw_s, maxomg_s, temp,  zbar, phi )
!call sub_tau( zbar,  phi, numw_s, maxomg_s, numw, temp,  tau )

!print*, 'self and sigma calculation done'
!open(21,file='test.dat')
!do n = 1,numw
!  write(21,*) n , tau(n)
!end do
!close(21)

! creating rho matrix (equivalent to meshgrid(rho) in matlab)
do i = 1,numw
   rho_mat1(i,:) = rho(:)
end do
do i = 1,numw
do j = 1,numw  
   if( i==j )  rho_mat1(i,j) = rho_mat1(i,j) + diff
end do
end do

do i = 1,numw
   rho_mat2(i,:) = rho(:)
end do
do i = 1,numw
do j = 1,numw
   if( i==j )  rho_mat2(i,j) = rho_mat2(i,j) - diff
end do
end do


kernel = 0.0d0
! constructing kernel(i,j)
do i = 1+rank,numw,numpc
 
    print*, '[ i, # of grid ] =',    i,   numw

    zbar = 0.0d0
     phi = 0.0d0
 
  call sub_self( rho_mat1(i,:), numw, numw_s, maxomg_s, temp,  zbar, phi )
  call sub_tau( zbar,  phi, numw_s, maxomg_s, numw, temp,  tau1 )
 
    zbar = 0.0d0
     phi = 0.0d0
 
  call sub_self( rho_mat2(i,:), numw, numw_s, maxomg_s, temp,  zbar, phi )
  call sub_tau( zbar,  phi, numw_s, maxomg_s, numw, temp,  tau2 )
   
    dtau = (tau1-tau2)/diff/dw  * 0.5d0  ! derivative from mid-point averaging
    kernel(:,i) = dtau

end do

call mpi_barrier( mpi_comm_world, ierr  )
call  mpi_allreduce( mpi_in_place, kernel, size(kernel), mpi_double_precision, &
        mpi_sum,     mpi_comm_world, ierr)
call mpi_barrier( mpi_comm_world, ierr  )


if( rank==0 )then
print*,  '   producing output file ... '
open(11,file='kernel-dw00X-tempXXK.dat')
do i = 1,numw
do j = 1,numw
    write(11,*) i, j, kernel(j,i)
end do
end do
close(11)
print*,  '   Done!!    '
end if



call  mpi_barrier(MPI_COMM_WORLD,ierr)
call mpi_finalize(ierr)
stop
end program main


