program spec
implicit none

!------------------------------------------------------------------------
! This will produce output files for spectral density 
! (Use these for mycond.f95 and cond1.f95)
! See Carbotte, PRB 71, 054506 (2005) Eq.17
!------------------------------------------------------------------------

double precision, parameter :: scal = 0.5d0   !***** change this ********
!  lam_mmp = -1.72*lam_phi_bos+1.86
!  lam_Z_bos = 0.5

double precision, parameter :: w0 = 0.041d0  !spin fluc energy
integer, parameter :: numv = 400   ! # of nnu (see numv 'model.f95)
double precision, parameter :: maxv = 0.4d0  ! max nu
double precision, parameter :: del = 0.005d0 !small parameter in delta fun

double precision :: dv, nu, lamzmm, lamzph, lampmm, lampph, pi, ztophi
double precision, dimension(1:numv) ::  Bzmm, Bpph, Bzph, Bpmm
integer :: i, j, l, m


open(11,file='a2f.OUT')

dv = maxv/dfloat(numv)
pi = 3.141592654d0

lamzmm = 0.0d0  ! z=Z(self energy), mm=MMPmodel 
lamzph = 0.0d0  ! ph=delta-function model (for some boson)
lampmm = 0.0d0  ! p=Phi
lampph = 0.0d0
do m = 1,numv       !  finding normalization constant
        nu = dfloat(m)*dv
        Bzmm(m) = nu/w0/( nu*nu + w0*w0 )
        Bzph(m) = del/pi/( del*del + (nu-w0)**2 )
        Bpmm(m) = nu/w0/( nu*nu + w0*w0 )
        Bpph(m) = del/pi/( del*del + (nu-w0)**2 )
        lamzmm = lamzmm + 2.0d0*Bzmm(m)*dv/nu
        lamzph = lamzph + 2.0d0*Bzph(m)*dv/nu
        lampmm = lampmm + 2.0d0*Bpmm(m)*dv/nu
        lampph = lampph + 2.0d0*Bpph(m)*dv/nu
end do

!----------  change here for lambda  ---------------------------------       
Bzmm = 1.0d0*Bzmm/lamzmm  ! These factors are the dimensionless constant
Bzph = 0.5d0*Bzph/lamzph  ! lamda. Feel free to change  
Bpmm = 1.0d0*Bpmm/lampmm
Bpph = scal*Bpph/lampph
!--------------------------------------------------------------------
lamzmm = 0.0d0
lamzph = 0.0d0
lampmm = 0.0d0
lampph = 0.0d0
do m = 1,numv
        nu = dfloat(m)*dv
        lamzmm = lamzmm + 2.0d0*Bzmm(m)*dv/nu
        lamzph = lamzph + 2.0d0*Bzph(m)*dv/nu
        lampmm = lampmm + 2.0d0*Bpmm(m)*dv/nu
        lampph = lampph + 2.0d0*Bpph(m)*dv/nu
        write(11,*) nu, Bzmm(m)+ Bzph(m), Bpmm(m)+Bpph(m)
end do

print*, 'lamzmm =', lamzmm, 'lamzph =', lamzph
print*, 'lampmm =', lampmm, 'lampph =', lampph
close(11)

end program spec












