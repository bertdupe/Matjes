module m_rw_MC
use m_input_types, only: MC_input
implicit none

private
public :: rw_MC

contains

!subroutine rw_MC(n_Tsteps,n_sizerelax,n_thousand,restart_MC_steps,Total_MC_Steps,T_relax,T_auto,cone,i_restart,ising,underrel,overrel,sphere,equi,print_relax,Cor_log)
subroutine rw_MC(inp_MC)
    use m_io_utils
    use m_io_files_utils
    type(MC_input),intent(out)  ::  inp_MC
    ! internal
    logical :: ising,underrel,overrel,sphere,equi
    logical :: methods(5)
    integer :: io_input,i
    
    methods=.False.
    
    io_input=open_file_read('input')
    
    call get_parameter(io_input,'input','n_Tsteps',inp_MC%n_Tsteps)
    call get_parameter(io_input,'input','n_sizerelax',inp_MC%n_sizerelax)
    call get_parameter(io_input,'input','n_relaxation',inp_MC%n_thousand)
    call get_parameter(io_input,'input','restart_MC_steps',inp_MC%restart_MC_steps)
    call get_parameter(io_input,'input','Total_MC_Steps',inp_MC%Total_MC_Steps)
    call get_parameter(io_input,'input','T_relax',inp_MC%T_relax)
    call get_parameter(io_input,'input','T_auto',inp_MC%T_auto)
    
    call get_parameter(io_input,'input','cone',inp_MC%cone)
    call get_parameter(io_input,'input','print_relax',inp_MC%print_relax)
    call get_parameter(io_input,'input','Cor_log',inp_MC%Cor_log)
    call get_parameter(io_input,'input','restart',inp_MC%i_restart)
    call get_parameter(io_input,'input','do_fluct',inp_MC%do_fluct)

    call get_parameter(io_input,'input','ising',inp_MC%ising)
    call get_parameter(io_input,'input','underrelaxation',inp_MC%underrelax)
    call get_parameter(io_input,'input','overrelaxation',inp_MC%overrelax)
    call get_parameter(io_input,'input','sphere_samp',inp_MC%sphere)
    call get_parameter(io_input,'input','equi_samp',inp_MC%equi)

    methods=[inp_mc%ising,inp_MC%underrelax,inp_MC%overrelax,inp_MC%sphere,inp_MC%equi]
    if(count(methods)>1)then
        write(*,*) "CANNOT SET MORE THAN ONE OF THE FOLLOWING SAMPLING METHODS AT ONCE:"
        write(*,*) "ising:",methods(1)
        write(*,*) "underrelaxation:",methods(2)
        write(*,*) "overrelaxation:",methods(3)
        write(*,*) "sphere_samp:",methods(4)
        write(*,*) "equi_samp:",methods(5)
        STOP
    elseif(count(methods)==0)then
        write(*,*) "NEED TO SET AT LEAST ONE OF THE MONTECARLO SAMPLING METHODS"
        STOP
    endif

    call close_file('input',io_input)
    
end subroutine


end module m_rw_MC
