module m_constants
implicit none

interface identity
  module procedure :: identity_3d,identity_Ndim
end interface identity



! Boltzmannfactor in eV/K
real(kind=8), Parameter :: k_B=0.000086173324d0
! magnetic constant
real(kind=8), Parameter :: mu_0=1.0d0
! bohr magneton, units eV/T
real(kind=8), Parameter :: mu_B=0.0000578838180661d0
! dielectric permeability of vacuum, units C**2.s**2/kg/m**3=C**2.m**2/J/m**3 <-missing s
! real(kind=8), Parameter :: epsilon_0=8.8541878128d-12
! dielectric permeability of vacuum, units C**2.nm**2/eV/nm**3
! dielectric permeability of vacuum, units C**2.m**2/J/m**3=1.0d18*C**2.nm**2/6.242d18/eV/1.0d27*nm**3
! dielectric permeability of vacuum, units C**2.nm**2/eV/nm**3
real(kind=8), Parameter :: epsilon_0=  1.418597282673602d-39 !1.41848571175905158603d-39
! speed of light, units nm/fs
real(kind=8), Parameter :: c=299.792458d0
!      real(kind=8), Parameter :: mu_B=1.0d0
! h in eV.s (wiki)
!      real(kind=8), parameter :: h=4.135667516d-15
! hbar in eV.fs
real(kind=8), parameter :: hbar=6.58211928d-1
!      real(kind=8), parameter :: hbar=1.0d0
!electron charge in Coulomb
real(8), parameter :: qel = 1.60217663d-19 !1.60217657d-19
real(8),parameter       :: pi= 3.14159265358979323846264338327950288d0

! basic unit convert
!!!!!!!!!!!!!!!!!!!!
! Energies
!!!!!!!!!!!!!!!!!!!!
real(8), parameter :: J_to_eV=6.242d18
real(8), parameter :: Ha_to_eV=27.2114d0
!!!!!!!!!!!!!!!!!!!!
! length
!!!!!!!!!!!!!!!!!!!!
real(8), parameter :: m_to_nm=1.0d9
real(8), parameter :: A_to_nm=0.1d0
real(8), parameter :: bohr_to_nm=0.0529177d0
!!!!!!!!!!!!!!!!!!!!
! masses
!!!!!!!!!!!!!!!!!!!!
real(8), parameter :: kg_to_amu=6.022d26
!!!!!!!!!!!!!!!!!!!!
! time
!!!!!!!!!!!!!!!!!!!!
real(8), parameter :: s_to_fs=1.0d15







contains

! function that spits out the unit matrix
function identity_3d(a)
    implicit none
    real(8), intent(in) :: a
    real(8), dimension(3,3) :: identity_3d
    real(8), dimension(3,3),parameter :: iden=reshape([1.0d0,0.0d0,0.0d0,0.0d0,1.0d0,0.0d0,0.0d0,0.0d0,1.0d0],[3,3])
    
    identity_3d = iden * a

end function identity_3d

function identity_Ndim(a,N)
    implicit none
    real(8), intent(in) :: a
    integer, intent(in) :: N
    real(8), dimension(N,N) :: identity_Ndim

    integer :: i

    identity_Ndim=0.0d0
    do i=1,N
       identity_Ndim(i,i)=a
    enddo

end function identity_Ndim

end module m_constants
