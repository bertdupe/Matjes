module m_shape_excitations

! variable that contains the shape of the fields
type shape_field
   real(kind=8), dimension(3) :: center=(/0.0d0,0.0d0,0.0d0/)
   real(kind=8) :: cutoff=1.0d0
   character(len=30) :: name='plane'
end type


abstract interface
    function shape_norm(R,R0,cutoff)
    real(kind=8), intent(in) :: R(3),R0(3),cutoff
    real(kind=8) :: shape_norm
  end function
end interface

type(shape_field), public, protected :: shape_excitation

private
public :: get_shape,shape_norm

contains

!
! subroutine that gets the shape of the field in the input file
!
subroutine get_shape(io,fname,name,norm)
implicit none
character(len=*), intent(in) :: fname,name
integer, intent(in) :: io
procedure(shape_norm), pointer, intent(out) :: norm
! internal

call read_shape(io,fname,'shape',shape_excitation)

select case (trim(shape_excitation%name))

  case('plane')
    norm => norm_0

  case('square')
    norm => norm_1

  case('cylinder')
    norm => norm_2

  case('gaussian')
    norm => norm_3

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

function norm_3(R,R0,cutoff)
use m_vector, only : norm
implicit none
real(kind=8), intent(in) :: R(3),R0(3),cutoff
real(kind=8) :: norm_3

norm_3=exp(-(norm(R-R0)/cutoff)**2)

end function norm_3



!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read the shape of the excitations for all the cycles
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
subroutine read_shape(io,fname,vname,field)
use m_io_utils, only : check_read
implicit none
type(shape_field), intent(inout) :: field
integer, intent(in) :: io
character(len=*), intent(in) :: vname,fname
! internal variable
integer :: fin,len_string,nread,check
character(len=100) :: str
character(len=100) :: dummy

nread=0
len_string=len(trim(adjustl(vname)))

field%name='plane'
field%center=(/0.0d0,0.0d0,0.0d0/)
field%cutoff=0.0d0

rewind(io)
do
   read (io,'(a)',iostat=fin) str
   if (fin /= 0) exit
   str= trim(adjustl(str))

   if (len_trim(str)==0) cycle
   if (str(1:1) == '#' ) cycle

!cccc We start to read the input
   if ( str(1:len_string) == trim(adjustl(vname)) ) then
      nread=nread+1
      backspace(io)
      read(io,*) dummy,field%name,field%center(1:3),field%cutoff
   endif

enddo

check=check_read(nread,vname,fname)

end subroutine read_shape

end module m_shape_excitations
