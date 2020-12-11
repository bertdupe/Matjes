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

type,extends(mpi_type)  ::   mpi_distv
    integer,allocatable :: cnt(:)       !how many are at each thread
    integer,allocatable :: displ(:)     !displacement for each thread
contains
    procedure   :: init => init_mpi_distv
end type

private :: init_mpi_distv
type(mpi_type)  ::  mpi_world

contains

subroutine init_mpi_distv(this,ref)
    class(mpi_distv),intent(inout)      :: this
    class(mpi_type),intent(in)          :: ref

    this%id   = ref%id   
    this%Np   = ref%Np   
    this%com  = ref%com  
    this%mas  = ref%mas  
    this%ismas= ref%ismas

    allocate(this%cnt(this%Np),source=0)
    allocate(this%displ(this%Np),source=0)
end subroutine



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
