module m_precision


!
! Arbitrary precision routine.
! It will put the component array of vector to 0 if they are smaller than erpsilon
!

interface truncate
  module procedure troncate_real
end interface

private
public :: truncate

contains

subroutine troncate_real(X,N)
implicit none
real(kind=8),intent(inout) :: X(:)
integer, intent(in) :: N
! internal
integer :: i

do i=1,N
   if (X(i).lt.1.0d-20) X(i)=0.0d0
enddo

end subroutine

end module m_precision
