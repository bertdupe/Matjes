module m_parameters_TB_IO_EF
!module of type which contains the Fermi energy calculation part of the tight-binding
implicit none
private
public parameters_TB_IO_EF

type parameters_TB_IO_EF
    real(8)     :: N_electrons=-1.0d0  !total number of electrons (negative value means use E_F) 
    real(8)     :: E_F_in=0.0d0        !fermi energy
    real(8)     :: kt=1.0d-3           !smearing factor in energy
contains
    procedure   :: read_file
    procedure   :: bcast => bcast_local
end type

contains
subroutine read_file(this,io,fname)
    use m_io_utils
    class(parameters_TB_IO_EF),intent(inout)    :: this
    integer,intent(in)                          :: io
    character(len=*), intent(in)                :: fname

    call get_parameter(io, fname, 'N_electrons', this%N_electrons)
    call get_parameter(io, fname, 'fermi_kt',    this%kt)
    call get_parameter(io, fname, 'TB_EF',       this%E_F_in)
end subroutine

subroutine bcast_local(this,comm)
    use mpi_basic
    use mpi_util
    class(parameters_TB_IO_EF),intent(inout)    ::  this
    type(mpi_type),intent(in)                   ::  comm
#ifdef CPP_MPI
    real(8)     :: arr(3)
    arr=[this%N_electrons,this%E_F_in,this%kt]
    Call bcast(arr,comm)
    this%N_electrons=arr(1)
    this%E_F_in     =arr(2)
    this%kt         =arr(3)
#else
    continue
#endif
end subroutine 

end module
