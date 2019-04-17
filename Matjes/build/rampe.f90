module m_rampe
use m_basic_types, only : excitations
use m_constants, only : k_b

private
public :: update_rampe

contains

subroutine update_rampe(time,kt,h_int,excite,counter,check)
implicit none
integer, intent(in) :: time
real(kind=8), intent(in) :: check(:)
type(excitations), intent(inout) :: excite(:)
integer, intent(inout) :: counter
real(kind=8), intent(inout) :: kt,h_int(:)
! internal


select case(excite(counter)%name)
  case('temp')

  if ((time.ge.excite(counter)%t_start).and.(time.le.excite(counter)%t_end)) then

     if (time.eq.excite(counter)%t_start) excite(counter)%start_value(1)=kt/k_b

     kT=kT+(excite(counter)%end_value(1)-excite(counter)%start_value(1))/real(excite(counter)%t_end-excite(counter)%t_start)*k_b
     if (abs(check(2)).gt.1.0d-8) write(6,'(a,2x,f16.6)') 'Final Temp', check(1)/check(2)/2.0d0/k_B
     write(6,'(a,f10.5)') 'T=',kT/k_B
  endif

  case('H_ext')

  if ((time.ge.excite(counter)%t_start).and.(time.le.excite(counter)%t_end)) then
     h_int=h_int+(excite(counter)%end_value-excite(counter)%start_value)/real(excite(counter)%t_end-excite(counter)%t_start+1)
     write(6,'(a,2x,a,2x,a,2x,f16.6)') 'external parameter ', excite(counter)%name, ' is now ',h_int
  endif

  case default
    write(6,'(a)') 'external parameter not implemented'
end select

end subroutine update_rampe

end module m_rampe
