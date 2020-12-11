module mpi_util
use mpi_basic
use m_MC_io, only: bcast_MC 
implicit none

interface scatter 
    module procedure scatter_int
end interface

interface bcast_alloc
    module procedure bcast_alloc_int2
end interface

interface bcast
    module procedure bcast_int
    module procedure bcast_int1
    module procedure bcast_real1
    module procedure bcast_MC
end interface

contains

subroutine scatter_int(arr,loc,com)
    use mpi_basic
    integer,intent(in)             :: arr(:)
    class(mpi_type),intent(in)     :: com
    integer,intent(out)            :: loc
#ifdef CPP_MPI    
    integer                        :: ierr

    Call MPI_SCATTER(arr(1),1,MPI_INTEGER,loc,1,MPI_INTEGER,com%mas,com%com,ierr)
#else
    continue
#endif
end subroutine

subroutine bcast_alloc_int2(arr,com)
    use mpi_basic
    integer,intent(inout),allocatable   :: arr(:,:)
    class(mpi_type),intent(in)          :: com
#ifdef CPP_MPI    
    integer     :: ierr
    integer     :: shp(2)

    if(com%ismas.and..not.allocated(arr)) ERROR STOP "FAILED TO BCAST SINCE INITIAL ARRAY IS NOT ALLOCATED"
    if(com%ismas) shp=shape(arr)
    Call bcast(shp,com)
    if(.not.allocated(arr)) allocate(arr(shp(1),shp(2)))

    Call MPI_BCAST(arr(1,1),size(arr),MPI_INTEGER,com%mas,com%com,ierr)
#else
    continue
#endif
end subroutine


subroutine bcast_int(val,com)
    integer,intent(inout)       :: val
    class(mpi_type),intent(in)  :: com
#ifdef CPP_MPI    
    integer     :: ierr

    Call MPI_BCAST(val,1,MPI_INTEGER,com%mas,com%com,ierr)
#else
    continue
#endif
end subroutine

subroutine bcast_int1(arr,com)
    integer,intent(inout)       :: arr(:)
    class(mpi_type),intent(in)  :: com
#ifdef CPP_MPI    
    integer     :: ierr

    Call MPI_BCAST(arr(1),size(arr),MPI_INTEGER,com%mas,com%com,ierr)
#else
    continue
#endif
end subroutine

subroutine bcast_real1(arr,com)
    real(8),intent(inout)       :: arr(:)
    class(mpi_type),intent(in)  :: com
#ifdef CPP_MPI    
    integer     :: ierr

    Call MPI_BCAST(arr(1),size(arr),MPI_DOUBLE_PRECISION,com%mas,com%com,ierr)
#else
    continue
#endif
end subroutine


end module
