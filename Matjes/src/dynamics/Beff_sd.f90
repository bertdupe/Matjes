module m_eval_Beff
use m_dipolar_field, only : i_dip,get_dipole_B
use m_internal_fields_commons, only : B_total

interface calculate_Beff
    module procedure normal
end interface calculate_Beff

private
public :: calculate_Beff
contains
! subroutine that calculates the field
! dE/DM
!
!
!
!--------------------------------------------------------------
! for normal
!
subroutine normal(B,iomp,spin,dim_mode)
use m_basic_types, only : vec_point
use m_modes_variables, only : point_shell_mode
use m_matrix, only : reduce
implicit none
! input variable
type(vec_point), intent(in) :: spin(:)
integer, intent(in) :: iomp,dim_mode
! output of the function
real(kind=8), intent(out) :: B(:)
! internals
integer :: N,i,j
real(kind=8) :: S_int(dim_mode)

!N=B_total%ncolumn
N=size(B_total%line(:,iomp))
B=0.0d0

do i=1,N

   j=B_total%line(i,iomp)

   call reduce(B_total%value(i,iomp),size(B_total%value(i,iomp)%order_op),S_int,spin(iomp)%w,spin(j)%w,dim_mode)

   B=B+S_int

enddo

if (i_dip) call get_dipole_B(B,iomp)

end subroutine normal

end module m_eval_Beff
