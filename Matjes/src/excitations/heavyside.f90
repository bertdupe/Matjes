module m_heavyside
use m_basic_types, only : excitations
use m_constants, only : k_b

private
public :: update_heavyside

contains

subroutine update_heavyside(time,kt,h_int,E_int,excite,counter,check)
implicit none
real(kind=8), intent(in) :: time
real(kind=8), intent(in) :: check(:)
type(excitations), intent(inout) :: excite(:)
integer, intent(inout) :: counter
real(kind=8), intent(inout) :: kt,h_int(:),E_int(:)
! internal

select case(excite(counter)%name)
  case('temp')

  if (time.eq.excite(counter)%t_start) then

     excite(counter)%start_value(1)=kt/k_b

     kT=excite(counter)%end_value(1)*k_b
     if (abs(check(2)).gt.1.0d-8) write(6,'(a,2x,f16.6)') 'initial Temp', check(1)/check(2)/2.0d0/k_B
     write(6,'(a,f10.5)') 'T=',kT/k_B
  endif

  if (time.eq.excite(counter)%t_end) then

     kT=excite(counter)%start_value(1)*k_b
     if (abs(check(2)).gt.1.0d-8) write(6,'(a,2x,f16.6)') 'final Temp', check(1)/check(2)/2.0d0/k_B
     write(6,'(a,f10.5)') 'T=',kT/k_B
  endif

  case('H_ext')

  if (time.eq.excite(counter)%t_start) then

     excite(counter)%start_value=h_int

     h_int=excite(counter)%end_value
     write(6,'(a,2x,a,2x,a,2x,f16.6)') 'external parameter ', excite(counter)%name, ' is now ',h_int
  endif

  if (time.eq.excite(counter)%t_end) then

     h_int=excite(counter)%start_value
     write(6,'(a,2x,a,2x,a,2x,f16.6)') 'external parameter ', excite(counter)%name, ' is now ',h_int
  endif

  case('E_ext')

  if (time.eq.excite(counter)%t_start) then

     excite(counter)%start_value=E_int

     E_int=excite(counter)%end_value
     write(6,'(a,2x,a,2x,a,2x,f16.6)') 'external parameter ', excite(counter)%name, ' is now ',E_int
  endif

  if (time.eq.excite(counter)%t_end) then

     E_int=excite(counter)%start_value
     write(6,'(a,2x,a,2x,a,2x,f16.6)') 'external parameter ', excite(counter)%name, ' is now ',E_int
  endif

  case default
    write(6,'(a)') 'external parameter not implemented'
end select

end subroutine update_heavyside

end module m_heavyside
