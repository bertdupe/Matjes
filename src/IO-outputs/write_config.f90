module m_write_config
!module with write_config, 
!which writes out current state of the different order parameters in sererate files
!, if the respective order-parameter is associated

!use m_io_utils
use m_io_files_utils, only: open_file_write,close_file
use m_io_utils, only: dump_config
use m_type_lattice,only: lattice
implicit none


interface write_config
 module procedure write_config_char_real
 module procedure write_config_char
 module procedure write_config_int
end interface

private
public write_config
contains
subroutine write_config_char(fname,lat)
    character(len=*),intent(in) ::  fname
    type(lattice),intent(in)    ::  lat

    integer     :: io
    character(len=100) :: filename

    if(associated(lat%M%all_modes))then
        write(filename,'(3A)')  'Mfield_',trim(adjustl(fname)),'.dat'
        io=open_file_write(filename)
        call dump_config(io,lat%M)
        call close_file(filename,io)
    endif

    if(associated(lat%E%all_modes))then
        write(filename,'(3A)')  'Efield_',trim(adjustl(fname)),'.dat'
        io=open_file_write(filename)
        call dump_config(io,lat%E)
        call close_file(filename,io)
    endif

    if(associated(lat%B%all_modes))then
        write(filename,'(3A)')  'Bfield_',trim(adjustl(fname)),'.dat'
        io=open_file_write(filename)
        call dump_config(io,lat%B)
        call close_file(filename,io)
    endif

    if(associated(lat%T%all_modes))then
        write(filename,'(3A)')  'Tfield_',trim(adjustl(fname)),'.dat'
        io=open_file_write(filename)
        call dump_config(io,lat%T)
        call close_file(filename,io)
    endif
end subroutine

subroutine write_config_int(i,lat)
    integer,intent(in)          ::  i
    type(lattice),intent(in)    ::  lat

    character(len=100) :: fname

    write(fname,*) i
    Call write_config_char(trim(adjustl(fname)),lat)
end subroutine

subroutine write_config_char_real(fname,r,lat,acc_op)
    character(len=*),intent(in) :: fname
    real(8),intent(in)          :: r
    type(lattice),intent(in)    :: lat
    integer,intent(in),optional :: acc_op(2)

    integer     ::  i,j
    integer,parameter   :: acc_par(2)=[6,4]
    integer ::  acc(2)

    character(len=100) :: comb_name
    character(len=100) :: name_format

    if(present(acc_op))then
        acc=acc_op
    else
        acc=acc_par
    endif
    write(name_format,*) '(A,A,I0.',acc(1),',A,I0.',acc(2),')'
    i=int(r)
    j=nint((r-real(int(r)))*10**(acc(2)))
    write(comb_name,name_format) Trim(adjustl(fname)),'_',i,'_',j
    Call write_config_char(trim(adjustl(comb_name)),lat)
end subroutine


end module
