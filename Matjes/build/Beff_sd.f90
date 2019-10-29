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
subroutine normal(B,spin,B_line)
use m_derived_types
implicit none
! input variable
type(point_shell_mode), intent(in) :: spin
type(point_shell_Operator), intent(in) :: B_line
! output of the function
real(kind=8), intent(out) :: B(:)
! internals
integer :: N,i

!N=B_total%ncolumn
N=size(B_line%shell)
B=0.0d0

do i=1,N

! the test takes more or less 10^-4s. Same time as the matmul
!   if (.not.associated(B_total%value(i,iomp)%Op_loc)) cycle

      B=B+matmul(B_line%shell(i)%Op_loc,spin%shell(i)%w)

!      if (iomp.eq.1) then
!      write(*,*) 'shell',i
!      write(*,*) spin%shell(i)%w
!      write(6,'(6(f12.6,2x))') (B_line%shell(i)%Op_loc(j,:),j=1,size(B_line%shell(i)%Op_loc,1))
!      write(*,*) matmul(B_line%shell(i)%Op_loc,spin%shell(i)%w)
!      write(*,*) ' '
!      endif

enddo

#ifdef CPP_DEBUG
      if (iomp.eq.1) write(*,*) B
      if (iomp.eq.1) pause
#endif
end subroutine normal

end module m_eval_Beff
