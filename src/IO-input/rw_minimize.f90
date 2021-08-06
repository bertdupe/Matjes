module m_rw_minimize
use m_input_types, only: min_input
implicit none
public

contains
subroutine rw_minimize(min_io,io_in,fname_in)
    use m_io_files_utils, only : open_file_read,close_file
    use m_io_utils,only: get_parameter
    type(min_input),intent(out)    :: min_io
    character(*),intent(in),optional    :: fname_in
    integer,intent(in),optional         :: io_in
    !integernal
    character(*),parameter              :: fname_default='input'
    character(:), allocatable           :: fname
    integer                             :: io_param

    if(present(fname_in))then
        fname=fname_in
    else
        fname=fname_default
    endif
    if(present(io_in))then
        io_param=io_in
    else
        io_param=open_file_read(fname)
    endif
    call get_parameter(io_param,fname,'steps',min_io%N_minimization)
    call get_parameter(io_param,fname,'masse',min_io%mass)
    call get_parameter(io_param,fname,'timestep',min_io%dt)
    call get_parameter(io_param,fname,'min_Efreq',min_io%Efreq)
    call get_parameter(io_param,fname,'convergence_criteria',min_io%conv_torque)

    if (min_io%mass.eq.0.0d0) then
        write(6,'(a)') 'The mass should be different from 0'
        stop
    endif

    if(.not.present(io_in)) call close_file(fname,io_param)
end subroutine
end module
