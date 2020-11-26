module m_precision

real(kind=8), parameter :: EPS=1.0d-50

!
! Arbitrary precision routine.
! It will put the component array of vector to 0 if they are smaller than erpsilon
!

interface truncate
  module procedure truncate_vecp_1D,truncate_lattice,truncate_real
  module procedure truncate_real_arr1d
  module procedure truncate_real_arr2d
  module procedure truncate_real_arr3d
  !needs  to be done for each rank as assumed shape does not work in interfaces...
end interface

private
public :: truncate

contains


subroutine truncate_lattice(lat,used)
use m_derived_types, only : lattice,number_different_order_parameters
implicit none
type(lattice),intent(inout) ::  lat
logical,intent(in)  :: used(number_different_order_parameters)

    if(used(1)) where(abs(lat%M%modes) < EPS) lat%M%modes=0.0d0
    if(used(2)) where(abs(lat%E%modes) < EPS) lat%E%modes=0.0d0
    if(used(3)) where(abs(lat%B%modes) < EPS) lat%B%modes=0.0d0
    if(used(4)) where(abs(lat%T%modes) < EPS) lat%T%modes=0.0d0
    if(used(5)) where(abs(lat%u%modes) < EPS) lat%u%modes=0.0d0

end subroutine

subroutine truncate_real_arr1d(arr,eps_in)
    real(8),intent(inout)       :: arr(:)
    real(8),optional,intent(in) :: eps_in

    real(8)     ::  cut
    if(present(eps_in))then
        cut=eps_in
    else
        cut=maxval(abs(arr))*EPS
    endif
    where(abs(arr)<cut) arr=0.0d0
end subroutine


subroutine truncate_real_arr2d(arr,eps_in)
    real(8),intent(inout)       :: arr(:,:)
    real(8),optional,intent(in) :: eps_in

    real(8)     ::  cut
    if(present(eps_in))then
        cut=eps_in
    else
        cut=maxval(abs(arr))*EPS
    endif
    where(abs(arr)<cut) arr=0.0d0
end subroutine

subroutine truncate_real_arr3d(arr,eps_in)
    real(8),intent(inout)       :: arr(:,:,:)
    real(8),optional,intent(in) :: eps_in

    real(8)     ::  cut
    if(present(eps_in))then
        cut=eps_in
    else
        cut=maxval(abs(arr))*EPS
    endif
    where(abs(arr)<cut) arr=0.0d0
end subroutine

subroutine truncate_real(X,N)
implicit none
real(kind=8),intent(inout) :: X(:)
integer, intent(in) :: N
! internal
integer :: i

do i=1,N
   if (abs(X(i)).lt.EPS) X(i)=0.0d0
enddo
end subroutine

subroutine truncate_vecp_1D(X,N)
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
