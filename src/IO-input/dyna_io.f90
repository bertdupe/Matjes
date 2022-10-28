module m_dyna_io
use mpi_basic
use m_basic_types, only : torque
private
public dyna_input

type dyna_input
    real(8) :: timestep=1.0d0       !full iteration time-step (fs)
    real(8) :: damping=1.0d0        !damping term magnitude (unitless)
    integer :: duration=1000        !number of full time-step iterated during simulation
    integer :: Efreq=100            !number of iterations after which state parameters are printed out (EM.dat)
    type(torque) :: SOT,STT
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
        ! SOT
        call get_parameter(io,fname,'SOT_current_dir',  this%SOT%current_dir)
        call get_parameter(io,fname,'SOT_polarization',  this%SOT%polarization)
        call get_parameter(io,fname,'SOT_FL_damp',  this%SOT%FL_damp)
        call get_parameter(io,fname,'SOT_DL_damp',  this%SOT%DL_damp)
        close(io)
    endif

    if(present(comm)) Call this%bcast(comm)
end subroutine
end module
