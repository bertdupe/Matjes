module m_io_wavefunc_evol
!module of type which contains the IO-parameters for the density of states (DOS) calculation
use m_dos_io
implicit none
private
public io_wavefunc_evol

type io_wavefunc_evol
    real(8) :: timestep=1.0d-3  !time step (fs)
    integer :: Efreq=100        !frequency of writing per time step
    integer :: duration=10000   !overall number of time steps
contains
    procedure :: read_file
    procedure :: bcast
end type


contains

subroutine read_file(this,fname,io)
    use m_io_read_util
    use m_io_utils
    use, intrinsic :: iso_fortran_env, only : output_unit
    class(io_wavefunc_evol),intent(inout)   :: this
    character(len=*), intent(in)            :: fname
    integer,intent(in)                      :: io

    call get_parameter(io,fname,'wave_tstep'    ,this%timestep)
    call get_parameter(io,fname,'wave_Efreq'    ,this%Efreq)
    call get_parameter(io,fname,'wave_duration' ,this%duration)
end subroutine 

subroutine bcast(this,comm)
    use mpi_util,only : mpi_type ,bcast_util=> bcast
    class(io_wavefunc_evol),intent(inout)   :: this
    type(mpi_type),intent(in)               :: comm

    Call bcast_util(this%timestep,comm) 
    Call bcast_util(this%Efreq   ,comm) 
    Call bcast_util(this%duration,comm) 
end subroutine
end module
