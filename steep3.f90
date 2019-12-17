program main
implicit none
include 'mpif.h'

!=======================================================================
! Maximum entropy using steepest descent method.  
! It will differentiate a2fz and a2fp if worked.
! compile w/ subroutine sub_self2 & sub_tau
! Functional derivative calculated from mid-point averaging.
!======================================================================

integer :: rank, numpc, ierr

double precision, parameter :: dw = 0.005d0  ! frequency grid sigma.OUT ****
double precision, parameter :: maxomg = 0.5d0 ! max omega from sigma.OUT *****
double precision, parameter :: temp = 10.0d0  ! temperature *****
double precision, parameter :: serror = 0.001d0  ! standard exp error
integer, parameter :: iter = 300  ! # of iteration

double precision, parameter :: maxomg_s = 1.0d0  ! max omega in self energy

double precision ::   del, diff, omg, free0, para, norm, pi  
double precision :: alpha1, alpha2, mval1, mval2, lam1, lam2, const
double precision, allocatable, dimension(:) :: tau_ex, dtau, a2fz, a2fp
double precision, allocatable, dimension(:) :: kai0, entropy0, tau0
double precision, allocatable, dimension(:) :: kai, entropy, df,  kai2
double complex, allocatable, dimension(:) :: zbar, phi
double precision, allocatable, dimension(:,:) :: a2fz_mat, a2fp_mat 
double precision, allocatable, dimension(:,:) :: a2fz_mat2, a2fp_mat2
double precision, allocatable, dimension(:) :: nu

integer :: i, j, n, m, w, num_it, numw, numw_s, ii

!---------------------------------------------------------------------------

call mpi_init(ierr)
call mpi_comm_rank(mpi_comm_world, rank, ierr)
call mpi_comm_size(mpi_comm_world, numpc, ierr)

del = 0.05d0     ! small # for updating a2f(w) 0.05   ******
diff = 0.02d0  !  small deviation in derivative 0.015 ******
pi = 3.141592654d0

alpha1 = 0.0d0  ! lagrange multiplier of shannon entropy  0.01
mval1 = 1.0d0  !  10.0d0**-7  model function in entropy treated constant
alpha2 = 0.0d0 
mval2 = mval1

numw   = int( maxomg/dw )
numw_s = int( maxomg_s/dw )

allocate(   dtau(1:numw), tau_ex(1:numw), a2fz(1:numw),  a2fp(1:numw)    )
allocate(   zbar(-numw_s:numw_s),  phi(-numw_s:numw_s)     )
allocate(   kai0(1:numw), entropy0(1:numw), tau0(1:numw)      )
allocate(   a2fz_mat(1:numw,1:numw), a2fp_mat(1:numw,1:numw)  )
allocate(   a2fz_mat2(1:numw,1:numw), a2fp_mat2(1:numw,1:numw)  )
allocate(   kai(1:numw), kai2(1:numw), df(1:2*numw)  )
allocate(   nu(1:numw)    )

open(11,file='tau_case1.dat')
!open(11,file='sigma.OUT')
do n = 1,numw
   nu(n) = dfloat(n)*dw
   ! read(11,*)  omg, para, para, para,  tau_ex(n)
     read(11,*)  omg, tau_ex(n)
end do


! initial guess of  a2f(w)
do n = 1,numw
   omg = dfloat(n)*dw
  !  a2fz(n) =  1.0d0 *  exp( -( omg/maxomg/0.3d0 ) ) /10.0d0          
     a2fz(n) = omg/0.04d0/( omg*omg + 0.04d0*0.04d0 )/3.0e3  ! mmp model
  !    a2fz(n) = 0.01d0*sin(2.0d0*pi/2.0d0/maxomg * omg )  
  ! if ( a2fz(n) < 0.0000001d0 )  a2fz(n) = 0.0000001d0
     a2fp(n) =   0.99d0 * a2fz(n)  
end do



!------------  iteration starting  --------------------------

if(rank==0)then
open(31,file='a2f-progress3.dat')
end if

call mpi_barrier( mpi_comm_world, ierr  )
num_it = 0
do ii = 1,iter

  if(rank==0)then
    print*, 'iteration = ', ii
  endif

  if( num_it > 100 ) then    ! update a2f, rho larger at beginning iteration
         del = del * 0.5d0    ! a2f(w)+del    
        ! diff = diff*1.5d0   ! rho+diff (for derivative)
        num_it = 1 
  end if


   101 continue
  call sub_self2( rank, a2fz, a2fp, numw, numw_s, maxomg_s, temp, zbar, phi )

   if( maxval( real(phi) ) < 0.001d0 ) then   ! check if a2fp is supercond
     a2fp = a2fp*1.1d0  
    if(rank==0) print*, ' a2fp(w) not superconducting '
     goto 101
   end if
   
    !  lam1 = 2.0d0*dw*sum( a2fz/nu )
    !  lam2 = 2.0d0*dw*sum( a2fp/nu )
     ! lamda_z needs to be bigger than lamba_phi 
    !   if( lam2 > lam1  )    a2fp = a2fp*lam1/lam2 *0.99d0
     ! ( notice 0.99 factor added since a2fz /= a2fp )

  call sub_tau(   zbar, phi, numw_s, maxomg_s, numw, temp, tau0 ) 
   

   !      kai0 = ( tau_ex - tau0 )*(  tau_ex - tau0 )/serror/serror
   !  entropy0 = alpha1*( a2fz - mval1 - a2fz*dlog( a2fz/mval1 ) )   &
    !             + alpha2*( a2fp - mval2 - a2fp*dlog( a2fp/mval2 ) ) 
    ! entropy0 = 0.0d0
     
   ! if( rank==0 ) then
   !  do n=1,numw
   !     if( isnan(tau0(n)) )      print*, ' tau0 = NaN '
   !     if( isnan(entropy0(n)) )  print*, ' entropy0 = NaN '  
   !  end do
   !  end if

    ! free0 = sum( kai0/2.0d0  - entropy0 )


! creating rho matrix (equivalent to meshgrid(rho) in matlab)
do i = 1,numw
   a2fz_mat(i,:) = a2fz(:)
   a2fp_mat(i,:) = a2fp(:)
end do
! adding diff to diagonal
do i = 1,numw
do j = 1,numw
   if( i==j )  a2fz_mat(i,j) = a2fz_mat(i,j) + diff
   if( i==j )  a2fp_mat(i,j) = a2fp_mat(i,j) + 1.0d0*diff  
end do
end do

do i = 1,numw
   a2fz_mat2(i,:) = a2fz(:)
   a2fp_mat2(i,:) = a2fp(:)
end do
! subtracting diff to diagonal
do i = 1,numw
do j = 1,numw
   if( i==j )  a2fz_mat2(i,j) = a2fz_mat2(i,j) - diff
   if( i==j )  a2fp_mat2(i,j) = a2fp_mat2(i,j) - 1.0d0*diff
end do
end do


df = 0.0d0
! Getting directional derivative of free energy
do n = 1+rank,numw,numpc

    call sub_self2( rank, a2fz_mat(n,:), a2fp, numw, numw_s, & 
                          maxomg_s, temp, zbar, phi )
    call sub_tau(   zbar, phi, numw_s, maxomg_s, numw, temp, dtau )
 
       kai = ( tau_ex - dtau )*(  tau_ex - dtau )/serror/serror
    
    call sub_self2( rank, a2fz_mat2(n,:), a2fp, numw, numw_s, &
                          maxomg_s, temp, zbar, phi )
    call sub_tau(   zbar, phi, numw_s, maxomg_s, numw, temp, dtau )

       kai2 = ( tau_ex - dtau )*(  tau_ex - dtau )/serror/serror

   if( rank ==0 ) then
   do w=1,numw
      if( isnan(dtau(w)) )     print*, ' dtau1 = NaN '
   !   if( isnan(entropy(w)) )  print*, 'entropy1=NaN' 
   end do
   end if

    df(n) = 0.5d0 * ( sum( kai/2.0d0 ) - sum( kai2/2.0d0 )  )

end do
do n = numw+1+rank,2*numw,numpc
   m = n-numw

   call sub_self2( rank, a2fz, a2fp_mat(m,:), numw, numw_s, &
                          maxomg_s, temp, zbar, phi )
   call sub_tau(  zbar, phi, numw_s, maxomg_s, numw, temp, dtau )
  
       kai = ( tau_ex - dtau )*(  tau_ex - dtau )/serror/serror

   call sub_self2( rank, a2fz, a2fp_mat2(m,:), numw, numw_s, &
                          maxomg_s, temp, zbar, phi )
   call sub_tau(  zbar, phi, numw_s, maxomg_s, numw, temp, dtau )

       kai2 = ( tau_ex - dtau )*(  tau_ex - dtau )/serror/serror

    do w=1,numw
        if( isnan(dtau(w)) )    print*, ' dtau2 = NaN '
    !    if( isnan(entropy(w)) ) print*, ' entropy2 = NaN ' 
    end do

    df(n) = 0.5d0 * ( sum( kai/2.0d0 ) - sum( kai2/2.0d0 )    )

end do

call  mpi_barrier( mpi_comm_world, ierr  )
call  mpi_allreduce( mpi_in_place, df, 2*numw, mpi_double_precision, &
        mpi_sum,     mpi_comm_world, ierr)
call mpi_barrier( mpi_comm_world, ierr  )


   if(rank==0)then
    do n=1,2*numw
     if( isnan(df(n)) )   print*, ' df(n) is NaN,   n = ' , n
    end do
   end if 

      ! norm = sqrt( sum( df*df )  )
    if( ii<50 ) then
         norm = sqrt( sum(  df(1:numw)*df(1:numw)   ) )
       a2fz = a2fz - del*df(1:numw)/norm
         norm = sqrt( sum( df(numw+1:2*numw)*df(numw+1:2*numw)  ) )
       a2fp = a2fp - 1.5d0 * del*df(numw+1:2*numw)/norm
    else
      norm = sqrt( sum( df*df )  )
       !  norm = sqrt( sum(  df(1:numw)*df(1:numw)   ) )
       a2fz = a2fz - del*df(1:numw)/norm
       !  norm = sqrt( sum( df(numw+1:2*numw)*df(numw+1:2*numw)  ) )
       a2fp = a2fp - 1.0d0 * del*df(numw+1:2*numw)/norm
    end if


  do n = 1,numw
   if ( a2fz(n) <=  0.0000001d0 )  a2fz(n) = 0.0000001d0
   if ( a2fp(n) <=  0.0000001d0 )  a2fp(n) = 0.0000001d0  
  
    if(rank==0)then
      if ( isnan( a2fz(n) ) )    print*,    ' a2fz = NaN ' 
      if ( isnan( a2fp(n) ) )    print*,    ' a2fp = NaN '   
       omg = dfloat(n)*dw
        write(31,*)     omg, a2fz(n), a2fp(n)
      if(n==numw)  write(31,*) ' '
    end if

  end do 

     num_it = num_it + 1

end do  

!-------------------  done  --------------------------------------


if(rank==0)then
open(21,file='a2f-steep-XXK3.dat')
open(22,file='sigma-model3.dat')
   call sub_self2( rank, a2fz, a2fp, numw, numw_s, &
                          maxomg_s, temp, zbar, phi )
   call sub_tau(  zbar, phi, numw_s, maxomg_s, numw, temp, tau0 )
do n = 1,numw
   omg = dfloat(n)*dw
  write(21,*)  omg, a2fz(n), a2fp(n)
  write(22,*)  omg, tau0(n), tau_ex(n) 
end do
close(21)
close(31)
close(22)
   print*, '  ****** Done ******   '
end if



call  mpi_barrier(MPI_COMM_WORLD,ierr)
call mpi_finalize(ierr)
stop
end program main


