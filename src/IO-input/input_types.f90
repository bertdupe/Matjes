module m_input_types
use m_constants, only : pi
!use m_MC_io, only : MC_input
implicit none
private :: pi
public

type extpar_input
    real(8) ::  H(3)=0.0d0      !constant magnetic field initializer
    real(8) ::  E(3)=0.0d0      !constant electric field initializer
    real(8) ::  T(2)=0.0d0      !initial (and final) magnetic temperature, only one temperature usually chooses T(1)
    logical :: enable_H=.false. !enable (switch on if not otherwise initialized) the usage of the magnetic field order parameter
    logical :: enable_E=.false. !enable (switch on if not otherwise initialized) the usage of the electric field order parameter
    logical :: enable_T=.false. !enable (switch on if not otherwise initialized) the usage of the temperature    order parameter
    logical :: enable_M=.false. !enable (switch on if not otherwise initialized) the usage of the magnetization  order parameter (nonsensical, as either there is a magnetic moment or not)
    logical :: enable_u=.false. !enable (switch on if not otherwise initialized) the usage of the displacement   order parameter
end type




end module
