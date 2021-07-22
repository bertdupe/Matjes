module m_parameters_TB_IO_OCC_mult
!module of type which contains the IO-parameters for the real-space occupation at several energies
implicit none
private
public parameters_TB_IO_OCC_MULT
type parameters_TB_IO_OCC_MULT
   !parameters for multiple occupation calculation at different energies
   real(8)                 ::  dE=-1.0                      !energy step size to plot
   real(8)                 ::  E_ext(2)=[0.0d0,0.0d0]      !minimal and aimed maximal energy considered
   real(8)                 ::  kt=-1.0d0                    !smearing factor of fermi and derivative of fermi function
contains
    procedure   :: read_file
    procedure   :: bcast => bcast_local
end type
contains
subroutine read_file(this,io,fname)
    use m_io_read_util
    use m_io_utils
    class(parameters_TB_IO_OCC_MULT),intent(inout)  :: this
    integer,intent(in)                              :: io
    character(len=*), intent(in)                    :: fname

    call get_parameter(io, fname, 'occ_mult_dE',        this%dE)
    call get_parameter(io, fname, 'occ_mult_E_ext',2,   this%E_ext)
    call get_parameter(io, fname, 'occ_mult_kt',        this%kt)
end subroutine

subroutine bcast_local(this,comm)
    use mpi_basic
    use mpi_util
    class(parameters_TB_IO_OCC_MULT),intent(inout)  :: this
    type(mpi_type),intent(in)                       :: comm
#ifdef CPP_MPI
    real(8)     :: arr(4)
    arr=[this%dE,this%E_ext(1),this%E_ext(2),this%kt]
    Call bcast(arr,comm)
    this%dE     =arr(1)
    this%E_ext  =arr(2:3)
    this%kt     =arr(4)
#else
    continue
#endif
end subroutine 

end module
