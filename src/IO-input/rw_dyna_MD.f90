subroutine rw_dyna_MD(timestep,Efreq,duration,file,damp_S,damp_F)
use m_constants
use m_derived_types
use m_io_files_utils
use m_io_utils
implicit none
real(kind=8), intent(out) :: timestep,damp_S,damp_F
integer, intent(out) :: duration,Efreq
character(len=100), intent(out) :: file
! internal
integer :: io
logical :: Ffield,i_Efield

Efreq=1
timestep=1.0d0
damp_S=0.0d0
damp_F=0.0d0
file='velocities.in'

io=open_file_read('input')

call get_parameter(io,'input','timestep',timestep)
call get_parameter(io,'input','Efreq',Efreq)
call get_parameter(io,'input','duration',duration)
call get_parameter(io,'input','velocities',file)
call get_parameter(io,'input','damp_Solid',damp_S)
call get_parameter(io,'input','damp_Fluid',damp_F)

Ffield=.false.
call get_parameter(io,'input','Ffield',Ffield)
i_Efield=.false.
call get_parameter(io,'input','Efield',i_Efield)

call close_file('input',io)

end subroutine rw_dyna_MD
