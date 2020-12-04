module mpi_basic
#ifdef CPP_MPI
use mpi
#endif
use,intrinsic :: iso_c_binding, only: c_bool
public
private :: set,c_bool
type        ::   mpi_type
    integer         :: id=0             !id
    integer         :: Np=1             !number processors
    integer         :: com=1            !communicator
    integer         :: mas=0            !master
    logical(c_bool) :: ismas=.true.     !is master
contains
    procedure   :: set
end type

type,extends(mpi_type)  :: mpi_type_distv
    !mpi type which contains how the values should be gatherv or scatterv'ed
    integer,allocatable :: displ(:)
    integer,allocatable :: cnt(:)
end type

type(mpi_type)  ::  mpi_world

contains


#ifdef CPP_MPI
subroutine set(this,com_in)
    class(mpi_type),intent(out)     :: this
    integer,intent(in)              :: com_in

    integer     ::  ierr

    this%com=com_in
    call MPI_COMM_RANK (com_in, this%id, ierr)
    if(ierr/=0) STOP "MPI_COMM_RANK failed"
    call MPI_COMM_SIZE (com_in, this%Np, ierr)
    if(ierr/=0) STOP "MPI_COMM_size failed"
    this%mas=0
    this%ismas=this%id==this%mas
end subroutine
#else
subroutine set(this) !no mpi
    class(mpi_type),intent(out)     :: this
    
    continue
    !without mpi just use default values

end subroutine
#endif
end module
