module m_eval_Beff
use m_dipolar_field, only : i_dip,get_dipole_B

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
subroutine normal(B,spin,B_line,iomp)
use m_derived_types
implicit none
! input variable
type(point_shell_mode), intent(in) :: spin
type(point_shell_Operator), intent(in) :: B_line
integer, intent(in) :: iomp
! output of the function
real(kind=8), intent(out) :: B(:)
! internals
integer :: N,i,j

!N=B_total%ncolumn
N=size(B_line%shell)
B=0.0d0

do i=1,N

! the test takes more or less 10^-4s. Same time as the matmul
!   if (.not.associated(B_total%value(i,iomp)%Op_loc)) cycle

      B=B+matmul(B_line%shell(i)%Op_loc,spin%shell(i)%w)

!      write(*,*) ''
!      write(*,*) spin%shell(i)%w
!      write(6,'(9(f12.6,2x))') (B_line%shell(i)%Op_loc(:,j),j=1,size(B_line%shell(i)%Op_loc,2))
!      write(*,*) matmul(B_line%shell(i)%Op_loc,spin%shell(i)%w)
!      pause

enddo

stop
if (i_dip) call get_dipole_B(B,iomp)

#ifdef CPP_DEBUG
      write(*,*) B,spin%shell(1)%w
      pause
#endif
end subroutine normal

end module m_eval_Beff
