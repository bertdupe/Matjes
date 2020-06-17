module m_precision

real(kind=8), parameter :: EPS=1.0d-50

!
! Arbitrary precision routine.
! It will put the component array of vector to 0 if they are smaller than erpsilon
!

interface truncate
  module procedure troncate_real,troncate_vecp_1D
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
   if (abs(X(i)).lt.EPS) X(i)=0.0d0
enddo

end subroutine


subroutine troncate_vecp_1D(X,N)
use m_basic_types, only : vec_point
implicit none
type(vec_point), intent(inout) :: X(:)
integer, intent(in) :: N
! internal
integer :: i,Size_M,j

Size_M=size(X)

do j=1,Size_M
  do i=1,N
   if (abs(X(j)%w(i)).lt.EPS) X(j)%w(i)=0.0d0
  enddo
enddo

end subroutine

end module m_precision