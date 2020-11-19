subroutine rw_dyna(timestep,damping,Efreq,duration)
use m_constants
use m_derived_types
use m_io_files_utils
use m_io_utils
implicit none
real(kind=8), intent(out) :: timestep,damping
integer, intent(out) :: duration,Efreq
! internal
integer :: io
logical :: Ffield,i_Efield,stmtemp

Efreq=1
timestep=1.0d0
damping=0.0d0

io=open_file_read('input')

call get_parameter(io,'input','timestep',timestep)
call get_parameter(io,'input','Efreq',Efreq)
call get_parameter(io,'input','duration',duration)
call get_parameter(io,'input','STMtemp',stmtemp)
call get_parameter(io,'input','damping',damping)

Ffield=.false.
call get_parameter(io,'input','Ffield',Ffield)
i_Efield=.false.
call get_parameter(io,'input','Efield',i_Efield)

call close_file('input',io)

end subroutine rw_dyna
