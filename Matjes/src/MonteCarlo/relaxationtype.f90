module m_relaxtyp
use m_derived_types
use m_eval_Beff

contains
! functions that relaxes the spins with respect of the dE/dM
! in one case, the spins are choosen in the direction of -dE/DM so energy diminishes
! in the second case, the spins are choosen in the direction of +dE/DM so energy increases
!
function underrelax(spin,B_line,iomp)
use m_vector, only : cross,norm
implicit none
! external variable
integer, intent(in) :: iomp
type(point_shell_mode), intent(in) :: spin(:)
type(point_shell_Operator), intent(in) :: B_line(:)
! value of the function
real(kind=8), dimension(3) :: underrelax
!internal variable
real(kind=8), dimension(3) ::S_int
real(kind=8) :: norm_local,dumy(3)


call calculate_Beff(S_int,spin(iomp),B_line(iomp),iomp)
norm_local=norm(S_int)
! Calculation of the new spin
!      norm=dsqrt((S_int(1)+spin_in(1))**2+(S_int(2)+spin_in(2))**2+(S_int(3)+spin_in(3))**2)
!      underrelax=(S_int+spin_in)/norm

if (norm_local.gt.1.0d-8) then
   dumy=cross(spin(iomp)%shell(1)%w,S_int,1,3)
   S_int=spin(iomp)%shell(1)%w-cross(spin(iomp)%shell(1)%w,dumy,1,3)
   norm_local=norm(S_int)
   underrelax=S_int/norm_local
else
   underrelax=spin(iomp)%shell(1)%w
endif

end function underrelax

function overrelax(spin,B_line,iomp)
use m_vector, only : cross,norm
implicit none
! external variable
integer, intent(in) :: iomp
type(point_shell_mode), intent(in) :: spin(:)
type(point_shell_Operator), intent(in) :: B_line(:)
! value of the function
real(kind=8), dimension(3) :: overrelax
!internal variable
real(kind=8), dimension(3) ::S_int
real(kind=8) :: norm_local,dumy(3)

call calculate_Beff(S_int,spin(iomp),B_line(iomp),iomp)
norm_local=norm(S_int)

if (norm_local.gt.1.0d-8) then
   dumy=cross(spin(iomp)%shell(1)%w,S_int,1,3)
   S_int=spin(iomp)%shell(1)%w-cross(spin(iomp)%shell(1)%w,dumy,1,3)
   norm_local=norm(S_int)
   overrelax=-S_int/norm_local
else
   overrelax=-spin(iomp)%shell(1)%w
endif

end function overrelax
end module
