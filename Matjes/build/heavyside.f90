module m_heavyside
use m_basic_types, only : excitations

private
public :: update_heavyside

contains

subroutine update_heavyside(time,field,excite,counter)
implicit none
real(kind=8), intent(in) :: time
type(excitations), intent(inout) :: excite(:)
integer, intent(inout) :: counter
real(kind=8), intent(inout) :: field(:)
! internal

if ((time.gt.excite(counter)%t_start).and.(time.lt.excite(counter)%t_end)) field=excite(counter)%start_value(:)

if (time.gt.excite(counter)%t_end) field=excite(counter)%end_value(:)

end subroutine update_heavyside

end module m_heavyside
