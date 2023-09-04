module m_probability
use m_constants, only : k_B
use,intrinsic :: iso_fortran_env, only : output_unit, error_unit
use m_io_utils
use m_io_files_utils
use m_proba_base
implicit none

private
public :: Boltzmann,Fermi_dirac,bose_einstein !, Bose_einstein, Fermi_dirac

contains

subroutine Boltzmann(this,delta_Energy,kT)
class(proba_data), intent(inout) :: this
real(8), intent(in)              :: delta_Energy(:,:)
real(8), intent(in)              :: kT

if (kT.lt.1.0d-8) then
    this%Pdistrib=1.0d0
  else
    this%Energy=sum(delta_Energy,2)
    this%Pdistrib=exp(-this%Energy/kT)
endif

end subroutine

subroutine Fermi_dirac(this,delta_Energy,kT)
class(proba_data), intent(inout) :: this
real(8), intent(in)              :: delta_Energy(:,:)
real(8), intent(in)              :: kT

if (kT.lt.1.0d-8) then
    this%Pdistrib=0.5d0
  else
    this%Energy=sum(delta_Energy,2)
    this%Pdistrib=1.0d0/(exp(-this%Energy/kT)+1.0d0)
endif

end subroutine

subroutine bose_einstein(this,delta_Energy,kT)
class(proba_data), intent(inout) :: this
real(8), intent(in)              :: delta_Energy(:,:)
real(8), intent(in)              :: kT

if (kT.lt.1.0d-8) then
    this%Pdistrib=100.0d0
  else
    this%Energy=sum(delta_Energy,2)
    this%Pdistrib=1.0d0/(exp(-this%Energy/kT)-1.0d0)
endif

end subroutine

end module
