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

subroutine set_MPI_type(blocks,bnd_real,bnd_cmplx,bnd_int,val_out)
    !subroutine which allows to commit derived type of continuously aranged real,complex, and integer values
    integer,intent(in),dimension(2) :: bnd_real,bnd_cmplx,bnd_int   !boundaries 
    integer,intent(in)              :: blocks(:) !how many enties(array size) each type entry has
    integer,intent(out)             :: val_out !newly commited MPI-type 
#ifdef CPP_MPI    
    !internal
    integer                         :: N_entry
    integer(kind=MPI_ADDRESS_KIND)  :: displ(size(blocks))
    integer                         :: types(size(blocks))
    integer(kind=MPI_ADDRESS_KIND)  :: LB,extend_real,extend_int,extend_cmplx
    integer(kind=MPI_ADDRESS_KIND)  :: extend
    integer                         :: ierr,i
    integer                         :: test_type

    N_entry=size(blocks)
    test_type=abs(MPI_DOUBLE_PRECISION)+abs(MPI_DOUBLE_COMPLEX)+abs(MPI_INT)
    types=test_type
    if(all(bnd_real >0)) types(bnd_real (1):bnd_real (2))=MPI_DOUBLE_PRECISION
    if(all(bnd_cmplx>0)) types(bnd_cmplx(1):bnd_cmplx(2))=MPI_DOUBLE_COMPLEX
    if(all(bnd_int  >0)) types(bnd_int  (1):bnd_int  (2))=MPI_INT
    if(any(types==test_type)) ERROR STOP "DID NOT FILL ALL TYPES WITH BOUNDS SETTING MPI_TYPE"
    CALL MPI_TYPE_GET_EXTENT(MPI_INT,lb,extend_int,ierr)
    CALL MPI_TYPE_GET_EXTENT(MPI_DOUBLE_PRECISION,lb,extend_real,ierr)
    CALL MPI_TYPE_GET_EXTENT(MPI_DOUBLE_COMPLEX,lb,extend_cmplx,ierr)
    displ=0
    do i=1,N_entry-1
        if(types(i)==MPI_DOUBLE_PRECISION)then
            extend=extend_real
        elseif(types(i)==MPI_INT)then
            extend=extend_int
        elseif(types(i)==MPI_DOUBLE_COMPLEX)then
            extend=extend_cmplx
        else
            STOP "UNEXPECTED types"
        endif
        displ(i+1)=displ(i)+blocks(i)*extend
    enddo

    CALL MPI_TYPE_CREATE_STRUCT(N_entry,blocks,displ,types,val_out,ierr)
    CALL MPI_TYPE_COMMIT(val_out,ierr)
#else
    continue
#endif
end subroutine


end module
