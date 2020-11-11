module m_rw_MC

private
public :: rw_MC

contains

subroutine rw_MC(n_Tsteps,n_sizerelax,n_thousand,restart_MC_steps,Total_MC_Steps,T_relax,T_auto,cone,i_restart,ising,underrel,overrel,sphere,equi,print_relax,Cor_log)
use m_constants, only : pi
use m_io_utils
use m_io_files_utils
implicit none
integer, intent(out) :: n_Tsteps,n_sizerelax,n_thousand,restart_MC_steps,T_relax,T_auto,Total_MC_Steps
real(kind=8), intent(out) :: cone
logical, intent(out) :: i_restart,ising,underrel,overrel,sphere,equi,print_relax,Cor_log
! internal
integer :: io_input

n_Tsteps=10
n_sizerelax=1
n_thousand=1000
cone=pi
restart_MC_steps=0
i_restart=.False.
ising=.False.
print_relax=.False.
equi=.False.
sphere=.True.
overrel=.False.
underrel=.False.
Cor_log=.False.
T_relax=1
T_auto=1
Total_MC_Steps=1000

io_input=open_file_read('input')

call get_parameter(io_input,'input','n_Tsteps',n_Tsteps)
call get_parameter(io_input,'input','n_sizerelax',n_sizerelax)
call get_parameter(io_input,'input','n_relaxation',n_thousand)
call get_parameter(io_input,'input','restart_MC_steps',restart_MC_steps)
call get_parameter(io_input,'input','Total_MC_Steps',Total_MC_Steps)
call get_parameter(io_input,'input','T_relax',T_relax)
call get_parameter(io_input,'input','T_auto',T_auto)

call get_parameter(io_input,'input','cone',cone)

call get_parameter(io_input,'input','restart',i_restart)
call get_parameter(io_input,'input','ising',ising)
call get_parameter(io_input,'input','underrelaxation',underrel)
call get_parameter(io_input,'input','overrelaxation',overrel)
call get_parameter(io_input,'input','sphere_samp',sphere)
call get_parameter(io_input,'input','equi_samp',equi)
call get_parameter(io_input,'input','print_relax',print_relax)
call get_parameter(io_input,'input','Cor_log',Cor_log)

call close_file('input',io_input)

end subroutine


end module m_rw_MC
