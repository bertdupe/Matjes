module m_hessian_io
use mpi_basic
private
public hessian_input

type hessian_input
    real(8) :: truc_test=1.0d0       !full iteration time-step (fs)
contains
    procedure   :: read_file =>read_hessian
    procedure   :: bcast
end type
contains

subroutine bcast(this,comm)
    use mpi_util,only : bcast_util=> bcast
    class(hessian_input),intent(inout) :: this
    type(mpi_type),intent(in)       :: comm

    Call bcast_util(this%truc_test,comm)
end subroutine

subroutine read_hessian(this,fname_in,comm)
    use, intrinsic  ::  ISO_FORTRAN_ENV, only: error_unit
    use m_io_utils
    class(hessian_input),intent(inout)      :: this
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
        ! list the different variables after that
        call get_parameter(io,fname,'test', this%truc_test)
        close(io)
    endif

    if(present(comm)) Call this%bcast(comm)
end subroutine
end module
