module m_relaxtyp
use m_derived_types, only : point_shell_Operator,lattice
use m_modes_variables, only : point_shell_mode
use m_basic_types, only : vec_point
use m_eval_Beff

contains
! functions that relaxes the spins with respect of the dE/dM
! in one case, the spins are choosen in the direction of -dE/DM so energy diminishes
! in the second case, the spins are choosen in the direction of +dE/DM so energy increases
!
function underrelax(iomp,lat)
use m_vector, only : cross,norm
implicit none
! external variable
integer, intent(in) :: iomp
type(lattice),intent(in)    :: lat
!type(vec_point), intent(in) :: spin(:)
! value of the function
real(kind=8), dimension(3) :: underrelax
!internal variable
real(kind=8), dimension(3) ::S_int
real(kind=8) :: norm_local,dumy(3)
type(vec_point), pointer :: spin(:)

spin=>lat%ordpar%all_l_modes
call calculate_Beff(S_int,iomp,spin,size(spin(iomp)%w))
norm_local=norm(S_int)
! Calculation of the new spin
!      norm=dsqrt((S_int(1)+spin_in(1))**2+(S_int(2)+spin_in(2))**2+(S_int(3)+spin_in(3))**2)
!      underrelax=(S_int+spin_in)/norm

if (norm_local.gt.1.0d-8) then
   dumy=cross(spin(iomp)%w,S_int,1,3)
   S_int=spin(iomp)%w-cross(spin(iomp)%w,dumy,1,3)
   norm_local=norm(S_int)
   underrelax=S_int/norm_local
else
   underrelax=spin(iomp)%w
endif
nullify(spin)

end function underrelax

function overrelax(iomp,lat)
use m_vector, only : cross,norm
implicit none
! external variable
integer, intent(in) :: iomp
type(lattice),intent(in)    :: lat
! value of the function
real(kind=8), dimension(3) :: overrelax
!internal variable
real(kind=8), dimension(3) ::S_int
real(kind=8) :: norm_local,dumy(3)
type(vec_point), pointer :: spin(:)

spin=>lat%ordpar%all_l_modes

call calculate_Beff(S_int,iomp,spin,size(spin(iomp)%w))
norm_local=norm(S_int)

if (norm_local.gt.1.0d-8) then
   dumy=cross(spin(iomp)%w,S_int,1,3)
   S_int=spin(iomp)%w-cross(spin(iomp)%w,dumy,1,3)
   norm_local=norm(S_int)
   overrelax=-S_int/norm_local
else
   overrelax=-spin(iomp)%w
endif

nullify(spin)

end function overrelax
end module
