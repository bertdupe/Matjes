module m_constants
! Boltzmannfactor in eV/K
real(kind=8), Parameter :: k_B=0.000086173324d0
! magnetic constant
real(kind=8), Parameter :: mu_0=1.0d0
! bohr magneton, units eV/T
real(kind=8), Parameter :: mu_B=0.0000578838180661d0
! dielectric permeability of vacuum, units fs**2/nm**2
real(kind=8), Parameter :: epsilon_0=0.00001112650056053618d0
! speed of light, units nm/fs
real(kind=8), Parameter :: c=299.792458d0
!      real(kind=8), Parameter :: mu_B=1.0d0
! h in eV.s (wiki)
!      real(kind=8), parameter :: h=4.135667516d-15
! hbar in eV.fs
real(kind=8), parameter :: hbar=6.58211928d-1
!      real(kind=8), parameter :: hbar=1.0d0
!electron charge in Coulomb
real(kind=8), parameter :: qel=1.60217657d-19
contains

! function that defines the pi constant and multiplies it
real(kind=8) function pi(a)
implicit none
real(kind=8), intent(in) :: a
!     Pi
real(kind=8), Parameter :: pi_value=acos(-1.0d0)

pi=a*pi_value
end function pi

! function that spits out the unit matrix
function identity(a)
implicit none
real(kind=8), intent(in) :: a
real(kind=8), dimension(3,3) :: identity

real(kind=8), dimension(1:9) :: iden=(/1.0d0,0.0d0,0.0d0,0.0d0,1.0d0,0.0d0,0.0d0,0.0d0,1.0d0/)

identity = RESHAPE( iden, (/ 3, 3 /) ) * a

end function identity

end module m_constants
