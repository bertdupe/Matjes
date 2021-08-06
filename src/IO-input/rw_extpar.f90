module m_rw_extpar
use m_input_types, only: extpar_input
implicit none
public
contains
subroutine rw_extpar(extpar_io,io_in,fname_in)
    use m_io_files_utils, only : open_file_read,close_file
    use m_io_utils,only: get_parameter
    type(extpar_input),intent(out)    ::  extpar_io
    integer,intent(in),optional         :: io_in
    character(*),intent(in),optional    :: fname_in
    !internal
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
    call get_parameter(io_param,fname,'H_ext',3,extpar_io%H)
    call get_parameter(io_param,fname,'E_ext',3,extpar_io%E)
    call get_parameter(io_param,fname,'Tini',extpar_io%T(1))
    call get_parameter(io_param,fname,'Tfin',extpar_io%T(2))
    call get_parameter(io_param,fname,'enable_H',extpar_io%enable_H)
    call get_parameter(io_param,fname,'enable_E',extpar_io%enable_E)
    call get_parameter(io_param,fname,'enable_T',extpar_io%enable_T)
    call get_parameter(io_param,fname,'enable_M',extpar_io%enable_M)
    call get_parameter(io_param,fname,'enable_u',extpar_io%enable_u)
    call get_parameter(io_param,fname,'enable_phonon',extpar_io%enable_u)
    if(.not.present(io_in)) call close_file(fname,io_param)
end subroutine
end module
