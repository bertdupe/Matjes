module m_FTeff_MD

use m_randist
use m_random_number_library
use m_constants, only : hbar

private
public :: Bolzmann
contains

subroutine Bolzmann(kt,damping,masses,FT)
implicit none
real(kind=8), intent(in) :: kt,damping,masses(:)
real(kind=8), intent(inout) :: FT(:)
! internal
integer :: i

do i=1,size(FT)
   FT(i)=sqrt(damping*masses(i))*randist(kt)
enddo

end subroutine Bolzmann

end module m_FTeff_MD
