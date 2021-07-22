module m_parameters_TB_IO_kint
!module of type which contains the integration kmesh part of the tight-binding
implicit none
private
public parameters_TB_IO_kint
type parameters_TB_IO_kint
    !parameters for getting a kmesh with energies within a certain energy window [Ecut(1),Ecut(2)]
    real(8)     :: Ecut(2)=0.d0
    integer     :: grid(3)=0
contains
    procedure   :: read_file
    procedure   :: bcast => bcast_local
end type
contains

subroutine read_file(this,io,fname)
    use m_io_read_util
    use m_io_utils
    class(parameters_TB_IO_kint),intent(inout)  :: this
    integer,intent(in)                          :: io
    character(len=*), intent(in)                :: fname
    
    call get_parameter(io,fname,'kint_grid',3,this%grid)
    call get_parameter(io,fname,'kint_Ecut',2,this%Ecut)
end subroutine

subroutine bcast_local(this,comm)
    use mpi_basic
    use mpi_util
    class(parameters_TB_IO_kint),intent(inout)  ::  this
    type(mpi_type),intent(in)                   ::  comm
#ifdef CPP_MPI
    Call bcast(this%Ecut ,comm)
    Call bcast(this%grid ,comm)
#else
    continue
#endif
end subroutine 

end module
