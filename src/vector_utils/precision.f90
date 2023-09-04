module m_precision

real(kind=8), private :: EPS=1.0d-50
!
! Arbitrary precision routine.
! It will put the component array of vector to 0 if they are smaller than erpsilon
!

interface truncate
  module procedure truncate_lattice,truncate_real
  module procedure truncate_real_arr1d
  module procedure truncate_real_arr2d
  module procedure truncate_real_arr3d
  !needs  to be done for each rank as assumed shape does not work in interfaces...
end interface

private
public :: truncate,set_EPS

contains


subroutine set_EPS()
use m_io_utils
use m_io_files_utils
implicit none

integer :: io_input

io_input=open_file_read('input')

! io variables
call get_parameter(io_input,'input','truncate_EPS',EPS)

call close_file('input',io_input)

end subroutine

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

end module m_precision
