module m_heavyside

private
public :: update_heavyside

contains

subroutine update_heavyside(time,field,t_start,t_end,start_value,end_value,counter)
implicit none
real(kind=8), intent(in) :: time,t_start,t_end
real(kind=8), intent(in) :: start_value(:),end_value(:)
integer, intent(inout) :: counter
real(kind=8), intent(inout) :: field(:)
! internal

if ((time.gt.t_start).and.(time.lt.t_end)) field=start_value(:)

if (time.gt.t_end) field=end_value(:)

end subroutine update_heavyside

end module m_heavyside
