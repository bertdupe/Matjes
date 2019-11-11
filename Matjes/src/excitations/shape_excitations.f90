module m_shape_excitations
use m_basic_types, only : shape_field


abstract interface
    function shape_norm(R,R0,cutoff)
    real(kind=8), intent(in) :: R(3),R0(3),cutoff
    real(kind=8) :: shape_norm
  end function
end interface

type(shape_field) :: shape_excitation

private
public :: get_shape

contains

!
! subroutine that gets the shape of the field in the input file
!
subroutine get_shape(io,fname,name,norm)
use m_io_utils, only : get_parameter
implicit none
character(len=*), intent(in) :: fname,name
procedure(shape_norm), pointer, intent(out) :: norm
integer, intent(in) :: io
! internal

call get_parameter(io,fname,'shape',shape_excitation)

select case (trim(shape_excitation%name))

  case('plane')
    norm => norm_0

  case('square')
    norm => norm_1

  case('cylinder')
    norm => norm_2

  case default
     stop 'The shape is not implemented'

end select

end subroutine get_shape




!
! The different norm functions
!
function norm_0(R,R0,cutoff)
implicit none
real(kind=8), intent(in) :: R(3),R0(3),cutoff
real(kind=8) :: norm_0

norm_0=1.0d0

end function norm_0

function norm_1(R,R0,cutoff)
implicit none
real(kind=8), intent(in) :: R(3),R0(3),cutoff
real(kind=8) :: norm_1

norm_1=0.0d0
if (all(abs(R-R0).lt.cutoff)) norm_1=1.0d0

end function norm_1

function norm_2(R,R0,cutoff)
use m_vector, only : norm
implicit none
real(kind=8), intent(in) :: R(3),R0(3),cutoff
real(kind=8) :: norm_2

norm_2=0.0d0
if (norm(R-R0).lt.cutoff) norm_2=1.0d0

end function norm_2

end module m_shape_excitations
