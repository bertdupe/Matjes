module m_rw_GNEB
use m_input_types, only: GNEB_input
implicit none
public
contains
subroutine rw_gneb(gneb_io,fname_in)
    use m_io_files_utils, only : open_file_read,close_file
    use m_io_utils,only: get_parameter
    type(GNEB_input),intent(out)        :: gneb_io
    character(*),intent(in),optional    :: fname_in
    !integernal
    character(*),parameter              :: fname_default='input'
    character(:), allocatable           :: fname
    integer                             :: io_param

    if(present(fname_in))then
        fname=fname_in
    else
        fname=fname_default
    endif
    io_param=open_file_read(fname)
    call get_parameter(io_param,fname,'momfile_i',gneb_io%momfile_i)
    call get_parameter(io_param,fname,'momfile_f',gneb_io%momfile_f)
    call get_parameter(io_param,fname,'restartfile_if',gneb_io%restartfile_if)
    call get_parameter(io_param,fname,'restartfile_path',gneb_io%restartfile_path)
    call get_parameter(io_param,fname,'spring',gneb_io%spring)
    call get_parameter(io_param,fname,'initpath',gneb_io%initpath)
    call get_parameter(io_param,fname,'amp_rnd',gneb_io%amp_rnd)
    call get_parameter(io_param,fname,'amp_rnd_path',gneb_io%amp_rnd_path)
    call get_parameter(io_param,fname,'min_itrmax',gneb_io%itrmax)
    call get_parameter(io_param,fname,'mintraj_step',gneb_io%mintraj_step)
    call get_parameter(io_param,fname,'min_ftol',gneb_io%minftol)
    call get_parameter(io_param,fname,'mep_itrmax',gneb_io%itrmax)
    call get_parameter(io_param,fname,'meptraj_step',gneb_io%every)
    call get_parameter(io_param,fname,'mep_ftol',gneb_io%mepftol)
    call get_parameter(io_param,fname,'mep_ftol_ci',gneb_io%mepftol_ci)
    call get_parameter(io_param,fname,'do_gneb_simple',gneb_io%do_gneb)
    call get_parameter(io_param,fname,'do_gneb_ci',gneb_io%do_gneb_ci)
    call get_parameter(io_param,fname,'do_norm_rx',gneb_io%do_norm_rx)
    call get_parameter(io_param,fname,'en_zero',gneb_io%en_zero)
    call get_parameter(io_param,fname,'vpo_dt',gneb_io%dt)
    call get_parameter(io_param,fname,'vpo_mass',gneb_io%mass)
    call get_parameter(io_param,fname,'sample_num',gneb_io%sample_num)
    call get_parameter(io_param,fname,'nim',gneb_io%nim)

    !gneb_io%init_path to support older input format
    select case(gneb_io%initpath)
    case(-1)
        continue !do nothing
    case(1)
        gneb_io%read_path=.False.      
        gneb_io%read_outer=.True.      
        gneb_io%min_type=0            
    case(2)
        gneb_io%read_path=.True.      
        gneb_io%read_outer=.True.      
        gneb_io%min_type=2            
    case(3)
        gneb_io%read_path=.True.      
        gneb_io%read_outer=.True.      
        gneb_io%min_type=1            
    case default
        gneb_io%read_path=.False.      
        gneb_io%read_outer=.True.      
        gneb_io%min_type=1            
    end select
    call get_parameter(io_param,fname,'gnebinit_read_path',gneb_io%read_path)
    call get_parameter(io_param,fname,'gnebinit_read_outer',gneb_io%read_outer)
    call get_parameter(io_param,fname,'gnebinit_min_type',gneb_io%min_type)

    call close_file(fname,io_param)
end subroutine
end module
