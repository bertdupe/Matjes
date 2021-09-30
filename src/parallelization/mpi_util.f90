module mpi_util
use mpi_basic
implicit none

interface scatter 
    module procedure scatter_int
end interface

interface bcast_alloc
    module procedure bcast_alloc_int
    module procedure bcast_alloc_int2
    module procedure bcast_alloc_real1
    module procedure bcast_alloc_real2
    module procedure bcast_alloc_cmplx2
    module procedure bcast_alloc_cmplx3
    module procedure bcast_alloc_character
    module procedure bcast_alloc_character0
    module procedure bcast_alloc_character_len
end interface

interface bcast
    module procedure bcast_int
    module procedure bcast_int1
    module procedure bcast_logical
    module procedure bcast_logical1
    module procedure bcast_real
    module procedure bcast_real1
    module procedure bcast_real2
    module procedure bcast_complex
    module procedure bcast_character
end interface

interface reduce_sum
    module procedure reduce_sum_int
    module procedure reduce_sum_int1
    module procedure reduce_sum_real1
    module procedure reduce_sum_real2
    module procedure reduce_sum_complex_arr
end interface 

interface reduce_lor
    module procedure reduce_lor_val
end interface 
contains


subroutine reduce_sum_int(val,com)
    use mpi_basic
    integer,intent(inout)          :: val   
    class(mpi_type),intent(in)     :: com
#ifdef CPP_MPI    
    integer                        :: ierr

    if(com%ismas)then
        Call MPI_Reduce(MPI_IN_PLACE, val, 1, MPI_INTEGER, MPI_SUM, com%mas, com%com, ierr)
    else
        Call MPI_Reduce(val,          val, 1, MPI_INTEGER, MPI_SUM, com%mas, com%com, ierr)
    endif
#else
    continue
#endif
end subroutine

subroutine reduce_sum_int1(arr,com)
    use mpi_basic
    integer,intent(inout)          :: arr(:)
    class(mpi_type),intent(in)     :: com
#ifdef CPP_MPI    
    integer                        :: ierr

    if(com%ismas)then
        Call MPI_Reduce(MPI_IN_PLACE, arr(1), size(arr), MPI_INTEGER, MPI_SUM, com%mas, com%com, ierr)
    else
        Call MPI_Reduce(arr(1),       arr(1), size(arr), MPI_INTEGER, MPI_SUM, com%mas, com%com, ierr)
    endif
#else
    continue
#endif
end subroutine


subroutine reduce_sum_real1(arr,com)
    use mpi_basic
    real(8),intent(inout)          :: arr(:)
    class(mpi_type),intent(in)     :: com
#ifdef CPP_MPI    
    integer                        :: ierr

    if(com%ismas)then
        Call MPI_Reduce(MPI_IN_PLACE, arr(1), size(arr), MPI_DOUBLE_PRECISION, MPI_SUM, com%mas, com%com, ierr)
    else
        Call MPI_Reduce(arr(1),       arr(1), size(arr), MPI_DOUBLE_PRECISION, MPI_SUM, com%mas, com%com, ierr)
    endif
#else
    continue
#endif
end subroutine

subroutine reduce_sum_real2(arr,com)
    use mpi_basic
    real(8),intent(inout)          :: arr(:,:)
    class(mpi_type),intent(in)     :: com
#ifdef CPP_MPI    
    integer                        :: ierr

    if(com%ismas)then
        Call MPI_Reduce(MPI_IN_PLACE, arr(1,1), size(arr), MPI_DOUBLE_PRECISION, MPI_SUM, com%mas, com%com, ierr)
    else
        Call MPI_Reduce(arr(1,1),     arr(1,1), size(arr), MPI_DOUBLE_PRECISION, MPI_SUM, com%mas, com%com, ierr)
    endif
#else
    continue
#endif
end subroutine

subroutine reduce_sum_complex_arr(arr,com)
    use mpi_basic
    complex(8),intent(inout)       :: arr(..)
    class(mpi_type),intent(in)     :: com
#ifdef CPP_MPI    
    integer                        :: ierr

    if(com%ismas)then
        Call MPI_Reduce(MPI_IN_PLACE, arr, size(arr), MPI_DOUBLE_COMPLEX, MPI_SUM, com%mas, com%com, ierr)
    else
        Call MPI_Reduce(arr,          arr, size(arr), MPI_DOUBLE_COMPLEX, MPI_SUM, com%mas, com%com, ierr)
    endif
#else
    continue
#endif
end subroutine



subroutine reduce_lor_val(val,com)
    use mpi_basic
    logical,intent(inout)          :: val
    class(mpi_type),intent(in)     :: com
#ifdef CPP_MPI    
    integer                        :: ierr

    if(com%ismas)then
        Call MPI_Reduce(MPI_IN_PLACE, val, 1, MPI_LOGICAL, MPI_LOR, com%mas, com%com, ierr)
    else
        Call MPI_Reduce(val,          val, 1, MPI_LOGICAL, MPI_LOR, com%mas, com%com, ierr)
    endif
#else
    continue
#endif
end subroutine


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

    shp=0
    if(com%ismas.and.allocated(arr)) shp=shape(arr)
    Call bcast(shp,com)
    if(all(shp>0))then
        if(.not.allocated(arr)) allocate(arr(shp(1),shp(2)))

        Call MPI_BCAST(arr(1,1),size(arr),MPI_INTEGER,com%mas,com%com,ierr)
    endif
#else
    continue
#endif
end subroutine

subroutine bcast_alloc_real1(arr,com)
    use mpi_basic
    real(8),intent(inout),allocatable   :: arr(:)
    class(mpi_type),intent(in)          :: com
#ifdef CPP_MPI    
    integer     :: ierr
    integer     :: shp

    shp=0
    if(com%ismas.and.allocated(arr)) shp=size(arr)
    Call bcast(shp,com)
    if(shp>0)then
        if(.not.allocated(arr)) allocate(arr(shp))

        Call MPI_BCAST(arr,size(arr),MPI_DOUBLE_PRECISION,com%mas,com%com,ierr)
    endif
#else
    continue
#endif
end subroutine

subroutine bcast_alloc_real2(arr,com)
    use mpi_basic
    real(8),intent(inout),allocatable   :: arr(:,:)
    class(mpi_type),intent(in)          :: com
#ifdef CPP_MPI    
    integer     :: ierr
    integer     :: shp(2)

    shp=0
    if(com%ismas.and.allocated(arr)) shp=shape(arr)
    Call bcast(shp,com)
    if(all(shp>0))then
        if(.not.allocated(arr)) allocate(arr(shp(1),shp(2)))

        Call MPI_BCAST(arr(1,1),size(arr),MPI_DOUBLE_PRECISION,com%mas,com%com,ierr)
    endif
#else
    continue
#endif
end subroutine

subroutine bcast_alloc_cmplx2(arr,com)
    use mpi_basic
    complex(8),intent(inout),allocatable   :: arr(:,:)
    class(mpi_type),intent(in)          :: com
#ifdef CPP_MPI    
    integer     :: ierr
    integer     :: shp(2)

    shp=0
    if(com%ismas.and.allocated(arr)) shp=shape(arr)
    Call bcast(shp,com)
    if(all(shp>0))then
        if(.not.allocated(arr)) allocate(arr(shp(1),shp(2)))

        Call MPI_BCAST(arr,size(arr),MPI_DOUBLE_COMPLEX,com%mas,com%com,ierr)
    endif
#else
    continue
#endif
end subroutine

subroutine bcast_alloc_cmplx3(arr,com)
    use mpi_basic
    complex(8),intent(inout),allocatable   :: arr(:,:,:)
    class(mpi_type),intent(in)          :: com
#ifdef CPP_MPI    
    integer     :: ierr
    integer     :: shp(3)

    if(com%ismas.and..not.allocated(arr)) ERROR STOP "FAILED TO BCAST SINCE INITIAL ARRAY IS NOT ALLOCATED"
    if(com%ismas) shp=shape(arr)
    Call bcast(shp,com)
    if(all(shp>0))then
        if(.not.allocated(arr)) allocate(arr(shp(1),shp(2),shp(3)))

        Call MPI_BCAST(arr(1,1,1),size(arr),MPI_DOUBLE_COMPLEX,com%mas,com%com,ierr)
    endif
#else
    continue
#endif
end subroutine

subroutine bcast_alloc_character_len(arr,len_in,com)
    use mpi_basic
    integer,intent(in)                              :: len_in
    character(len=len_in),intent(inout),allocatable :: arr(:)
    class(mpi_type),intent(in)                      :: com
#ifdef CPP_MPI    
    integer     :: ierr
    integer     :: shp

    shp=0
    if(com%ismas.and.allocated(arr)) shp=size(arr)

    Call bcast(shp,com)
    if(shp>0)then
        if(.not.allocated(arr)) allocate(arr(shp))

        Call MPI_BCAST(arr,size(arr)*len(arr),MPI_CHARACTER,com%mas,com%com,ierr)
    endif
#else
    continue
#endif
end subroutine

subroutine bcast_alloc_int(arr,com)
    use mpi_basic
    integer,intent(inout),allocatable   :: arr(:)
    class(mpi_type),intent(in)          :: com
#ifdef CPP_MPI    
    integer     :: ierr
    integer     :: shp

    shp=0
    if(com%ismas.and.allocated(arr)) shp=size(arr)
    Call bcast(shp,com)
    if(shp>0)then
        if(.not.allocated(arr)) allocate(arr(shp))

        Call MPI_BCAST(arr(1),shp,MPI_INTEGER,com%mas,com%com,ierr)
    endif
#else
    continue
#endif
end subroutine


subroutine bcast_alloc_character0(arr,com)
    use mpi_basic
    character(len=:),intent(inout),allocatable      :: arr
    class(mpi_type),intent(in)                      :: com
#ifdef CPP_MPI    
    integer     :: ierr
    integer     :: shp

    shp=0
    if(com%ismas)  shp=len(arr)

    Call bcast(shp,com)
    if(shp>0)then
        if(.not.allocated(arr)) allocate(character(len=shp) :: arr)

        Call MPI_BCAST(arr,len(arr),MPI_CHARACTER,com%mas,com%com,ierr)
    endif
#else
    continue
#endif
end subroutine

subroutine bcast_alloc_character(arr,com)
    use mpi_basic
    character(len=:),intent(inout),allocatable      :: arr(:)
    class(mpi_type),intent(in)                      :: com
#ifdef CPP_MPI    
    integer     :: ierr
    integer     :: shp(2)

    shp=0
    if(com%ismas.and.allocated(arr))then
        shp=[len(arr),size(arr)]
        write(*,*) shp, "DASDA"
        STOP "FASFAS"
    endif

    Call bcast(shp,com)
    if(all(shp>0))then
        if(.not.allocated(arr)) allocate(character(len=shp(1)) :: arr(shp(2)))

        Call MPI_BCAST(arr,size(arr)*len(arr),MPI_CHARACTER,com%mas,com%com,ierr)
    endif
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

subroutine bcast_logical(val,com)
    logical,intent(inout)       :: val
    class(mpi_type),intent(in)  :: com
#ifdef CPP_MPI    
    integer     :: ierr

    Call MPI_BCAST(val,1,MPI_LOGICAL,com%mas,com%com,ierr)
#else
    continue
#endif
end subroutine


subroutine bcast_logical1(val,com)
    logical,intent(inout)       :: val(:)
    class(mpi_type),intent(in)  :: com
#ifdef CPP_MPI    
    integer     :: ierr

    Call MPI_BCAST(val(1),size(val),MPI_LOGICAL,com%mas,com%com,ierr)
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


subroutine bcast_real(val,com)
    real(8),intent(inout)       :: val
    class(mpi_type),intent(in)  :: com
#ifdef CPP_MPI    
    integer     :: ierr

    Call MPI_BCAST(val,1,MPI_DOUBLE_PRECISION,com%mas,com%com,ierr)
#else
    continue
#endif
end subroutine

subroutine bcast_character(val,com)
    character(len=*),intent(inout)  :: val
    class(mpi_type),intent(in)      :: com
#ifdef CPP_MPI    
    integer     :: length
    integer     :: ierr

    length=len(val)
    Call MPI_Allreduce( MPI_IN_PLACE, length, 1, MPI_INTEGER,MPI_MAX, com%com,ierr)
    if(length/=len(val)) ERROR STOP "Trying to broadcast characters of different length"
    Call MPI_Bcast(val, length, MPI_CHARACTER, com%mas, com%com,ierr)
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

subroutine bcast_real2(arr,com)
    real(8),intent(inout)       :: arr(:,:)
    class(mpi_type),intent(in)  :: com
#ifdef CPP_MPI    
    integer     :: ierr

    Call MPI_BCAST(arr(1,1),size(arr),MPI_DOUBLE_PRECISION,com%mas,com%com,ierr)
#else
    continue
#endif
end subroutine


subroutine bcast_complex(val,com)
    complex(8),intent(inout)    :: val
    class(mpi_type),intent(in)  :: com
#ifdef CPP_MPI    
    integer     :: ierr

    Call MPI_BCAST(val,1,MPI_DOUBLE_COMPLEX,com%mas,com%com,ierr)
#else
    continue
#endif
end subroutine


end module
