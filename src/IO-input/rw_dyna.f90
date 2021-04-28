module m_dyna_io
use mpi_basic
private
public dyna_input
public rw_dyna

type dyna_input
    real(8) :: timestep=1.0d0
    real(8) :: damping=1.0d0
    integer :: Efreq=100
    integer :: duration=1000
contains
    procedure   :: read_file =>read_dyna
    procedure   :: bcast
end type
contains

subroutine bcast(this,comm)
    use mpi_util,only : bcast_util=> bcast
    class(dyna_input),intent(inout) :: this
    type(mpi_type),intent(in)       :: comm

    Call bcast_util(this%timestep,comm) 
    Call bcast_util(this%damping ,comm) 
    Call bcast_util(this%Efreq   ,comm) 
    Call bcast_util(this%duration,comm) 
end subroutine

subroutine read_dyna(this,fname_in,comm)
    use, intrinsic  ::  ISO_FORTRAN_ENV, only: error_unit
    use m_io_utils
    class(dyna_input),intent(inout)         :: this
    character(len=*),intent(in),optional    :: fname_in  
    type(mpi_type),intent(in),optional      :: comm
    character(len=*),parameter              :: fname_default="input"

    character(len=:),allocatable            :: fname
    integer                 :: io
    logical                 :: fexist, ismas

    ismas=.true.
    if(present(comm)) ismas=comm%ismas
    if(ismas)then
        if(present(fname_in))then
            fname=fname_in
        else
            fname=fname_default
        endif
        inquire(file=fname,exist=fexist)
        if(.not.fexist)then
            write(error_unit,'(3/,A,/,2A)') "File required to read dynamics input not found","File name: ",fname
            Error STOP 
        endif
        open(newunit=io,file=fname)
        call get_parameter(io,fname,'timestep', this%timestep)
        call get_parameter(io,fname,'Efreq',    this%Efreq)
        call get_parameter(io,fname,'duration', this%duration)
        call get_parameter(io,fname,'damping',  this%damping)
        close(io)
    endif

    if(present(comm)) Call this%bcast(comm)
end subroutine

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

end module
