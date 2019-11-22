module m_rampe
use m_basic_types, only : excitations
use m_convert

private
public :: update_rampe

contains

subroutine update_rampe(time,field,excite,counter)
implicit none
real(kind=8), intent(in) :: time
type(excitations), intent(inout) :: excite(:)
integer, intent(inout) :: counter
real(kind=8), intent(inout) :: field(:)
! internal
character(len=30) :: form
integer :: n

n=size(field)

if ((time.ge.excite(counter)%t_start).and.(time.le.excite(counter)%t_end)) then

   field=field+(excite(counter)%end_value(1)-excite(counter)%start_value(1))/real(excite(counter)%t_end-excite(counter)%t_start)

   form=convert('(a,',n,'f14.6)')
   write(6,form) 'field value ',field
endif

end subroutine update_rampe

end module m_rampe
