module m_eval_Beff

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
subroutine normal(iomp,B,spin,h_int,B_line)
use m_internal_fields_commons, only : B_total
use m_derived_types
use m_external_fields, only : ext_field
implicit none
! input variable
integer, intent(in) :: iomp
type(point_shell_mode), intent(in) :: spin
real(kind=8), intent(in) :: h_int(3)
type(point_shell_Operator), intent(in) :: B_line
! output of the function
real(kind=8), intent(out) :: B(:)
! internals
real(kind=8) :: mu_s
logical :: i_dip
integer :: N,i

!N=B_total%ncolumn
N=size(B_line%shell)
B=0.0d0
if (size(ext_field(iomp)%w).ne.0) B=ext_field(iomp)%w

do i=1,N

! the test takes more or less 10^-4s. Same time as the matmul
!   if (.not.associated(B_total%value(i,iomp)%Op_loc)) cycle

      B=B+matmul(B_line%shell(i)%Op_loc,spin%shell(i)%w)

enddo

#ifdef CPP_DEBUG
      if (iomp.eq.1) write(*,*) B
      if (iomp.eq.1) pause
#endif
end subroutine normal

end module m_eval_Beff
