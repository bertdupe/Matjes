module m_rampe
use m_convert

private
public :: update_rampe

contains

subroutine update_rampe(time,field,t_start,t_end,start_value,end_value,counter)
implicit none
real(kind=8), intent(in) :: time,t_start,t_end
real(kind=8), intent(in) :: end_value(:),start_value(:)
integer, intent(inout) :: counter
real(kind=8), intent(inout) :: field(:)
! internal
character(len=30) :: form
integer :: n

n=size(field)

if ((time.ge.t_start).and.(time.le.t_end)) then

   field=field+(end_value-start_value)/real(t_end-t_start,8)

   form=convert('(a,',n,'f14.6)')
   write(6,form) 'field value ',field
endif

end subroutine update_rampe

end module m_rampe
