module m_mode_construction_rankN_sparse_col
use m_mode_construction
use m_derived_types, only : lattice,number_different_order_parameters
use m_coo_mat
implicit none
private
public F_mode_rankN_sparse_col

type col_mat
    integer :: dim_mat(2)=0
    integer :: nnz=0
    integer,allocatable :: col(:)
end type

type, extends(F_mode) :: F_mode_rankN_sparse_col 
    integer,allocatable         :: order(:)
    integer                     :: N_mode=-1
    type(col_mat),allocatable   :: dat(:)
    contains
    !necessary routines as defined by class
    procedure   :: get_mode   !subroutine which returns the mode 
    procedure   :: get_mode_exc
    procedure   :: mode_reduce  

    procedure   :: copy
    procedure   :: bcast
    procedure   :: destroy
    procedure   :: is_same
    !local construction routine
    procedure   :: init_order
end type

contains

subroutine get_mode_exc(this,lat,op_exc,vec)
    use, intrinsic :: iso_fortran_env, only : error_unit
    class(F_mode_rankN_sparse_col),intent(in)   :: this
    type(lattice),intent(in)                    :: lat       !lattice type which knows about all states
    integer,intent(in)                          :: op_exc !of which operator the first entry is kept
    real(8),intent(inout)                       :: vec(:)

    logical         :: exclude(size(this%order))
    integer         :: N
    real(8)         :: tmp_internal(this%N_mode,size(this%dat))
    real(8),pointer :: mode_base(:)
    integer         :: i,j

    if(size(vec)/=this%N_mode) STOP "mode exc call has wrong size for vector"
    exclude=.false.
    i=findloc(this%order,op_exc,dim=1)
    if(i<1.or.i>size(this%order))then
        write(error_unit,'(//A,I6)') "Tried to get mode excluding order no.:", op_exc
        write(error_unit,*) "But the mode only contains the order:", this%order
        ERROR STOP "This makes no sense and should probably prevented earlier in the code"
    endif
    exclude(i)=.true.

    tmp_internal=1.d0
    do i=1,size(this%dat)
        if(.not.exclude(i))then
            Call lat%set_order_point(this%order(i),mode_base)
            tmp_internal(:,i)=mode_base(this%dat(i)%col)
        endif
    enddo
    vec=product(tmp_internal,dim=2)
    nullify(mode_base)
end subroutine

subroutine mode_reduce(this,lat,vec_in,op_keep,vec_out)
    use, intrinsic :: iso_fortran_env, only : error_unit
    class(F_mode_rankN_sparse_col),intent(in)  :: this
    real(8),intent(in)                          :: vec_in(:)
    type(lattice),intent(in)                    :: lat       !lattice type which knows about all states
    integer,intent(in)                          :: op_keep   !of which operator the first entry is kept
    real(8),intent(out)                         :: vec_out(lat%dim_modes(op_keep)*lat%Ncell)

    integer     ::  i_order,i

    i_order=findloc(this%order,op_keep,dim=1)
    vec_out=0.0d0
    do i=1,this%dat(i_order)%nnz
        vec_out(this%dat(i_order)%col(i))=vec_out(this%dat(i_order)%col(i))+vec_in(i)
    enddo
end subroutine


subroutine get_mode(this,lat,mode,tmp)
    class(F_mode_rankN_sparse_col),intent(in)   :: this
    type(lattice),intent(in)                    :: lat       !lattice type which knows about all states
    real(8),intent(out),pointer                 :: mode(:)   !pointer to required mode
    real(8),allocatable,target,intent(inout)    :: tmp(:)

    integer         :: N
    real(8)         :: tmp_internal(this%N_mode,size(this%dat))
    real(8),pointer :: mode_base(:)

    integer         :: i,j

    allocate(tmp(this%N_mode),source=0.0d0)
    mode=>tmp

    tmp_internal=0.d0
    do i=1,size(this%dat)
        Call lat%set_order_point(this%order(i),mode_base)
        tmp_internal(:,i)=mode_base(this%dat(i)%col)
    enddo
    mode=product(tmp_internal,dim=2)
    nullify(mode_base)
end subroutine

function is_same(this,comp)result(same)
    class(F_mode_rankN_sparse_col),intent(in)       :: this
    class(F_mode),intent(in)                   :: comp
    logical                                    :: same

    ERROR STOP "IMPLEMENT"
    !same=.false.
    !select type(comp) 
    !type is(F_mode_rankN_sparse_col)
    !    same=all(this%order==comp%order)
    !end select
end function

subroutine destroy(this)
    class(F_mode_rankN_sparse_col),intent(inout) ::  this
    integer ::  i
    this%N_mode=-1
    do i=1,size(this%dat)
        deallocate(this%dat(i)%col)
    enddo
    deallocate(this%dat)
    deallocate(this%order)
end subroutine

subroutine copy(this,F_out)
    class(F_mode_rankN_sparse_col),intent(in)    :: this
    class(F_mode),allocatable,intent(inout) :: F_out

    integer ::  i

    if(.not.allocated(F_out)) allocate(F_mode_rankN_sparse_col::F_out)
    select type(F_out)
    class is(F_mode_rankN_sparse_col)
        if(.not.allocated(F_out%order)) allocate(F_out%order(size(this%order)))
        if(size(F_out%order)/=size(this%order)) ERROR STOP "CANNOT COPY ORDER AS RANKS DIFFER"
        F_out%order=this%order
        F_out%N_mode=this%N_mode
        allocate(F_out%dat(size(this%dat)))
        do i=1,size(this%dat)
            F_out%dat(i)%dim_mat=this%dat(i)%dim_mat
            F_out%dat(i)%nnz=this%dat(i)%nnz
            allocate(F_out%dat(i)%col,source=this%dat(i)%col)
        enddo
    class default
        ERROR STOP "FAILED TO COPY F_mode_rankN_sparse_col mode to F_out"
    end select
end subroutine

subroutine bcast(this,comm)
    use mpi_basic                
    class(F_mode_rankN_sparse_col),intent(inout) ::  this        !this might fail if the server threads non-allocated class(F_mode), TAKE CARE OF THIS IN HAM_BASE
    type(mpi_type),intent(in)               ::  comm
#ifdef CPP_MPI
    integer     :: ierr
    integer     ::  N
  
    !THIS MIGHT BE INSUFFICIENT, MAYBE ONE HAS TO CHECK IF THE F_MODE IS ALREADY ALLOCATED TO THE F_mode_rankN_sparse_col type
    STOP "CHECK IF THIS WORKS WITHOUT PREVIOUS ALLOCATION./type stuff/, on non-master threads"
    if(comm%ismas)then
        if(.not.allocated(this%order)) ERROR STOP "CANNOT BCAST SINCE MASTER ORDER IS NOT ALLOCATED"
        N=size(this%order)
    endif
    Call MPI_Bcast(N,1, MPI_INTEGER, comm%mas, comm%com,ierr)
    if(.not.allocated(this%order)) allocate(this%order(N))
    Call MPI_Bcast(this%order,N, MPI_INTEGER, comm%mas, comm%com,ierr)
    Call MPI_Bcast(this%N_mode,1, MPI_INTEGER, comm%mas, comm%com,ierr)

    ERROR STOP "ALSO TOTALLY UNTESTED"
    !bcast mat
    do i=1,N
        Call MPI_Bcast(this%dat(i)%nnz    ,1, MPI_INTEGER, comm%mas, comm%com,ierr)
        Call MPI_Bcast(this%dat(i)%dim_mat,2, MPI_INTEGER, comm%mas, comm%com,ierr)
        if(.not.comm%ismas)then
            allocate(this%dat(i)%col(this%mat(i)%nnz))
        endif
        Call MPI_Bcast(this%dat(i)%col,this%dat(i)%nzz, MPI_INTEGER         , comm%mas, comm%com,ierr)
    enddo
#else
    continue
#endif
end subroutine 

subroutine init_order(this,lat,abbrev_in,mat)
    use m_derived_types, only: op_abbrev_to_int
    class(F_mode_rankN_sparse_col),intent(inout) :: this
    type(lattice),intent(in)                :: lat       !lattice type which knows about all states
    character(len=*), intent(in)            :: abbrev_in !considered order abbreviations
    type(coo_mat),intent(inout)             :: mat(:)    !input matrices, destroyed when returned
    integer     :: order(len(abbrev_in))
    integer     :: i     

    order=op_abbrev_to_int(abbrev_in)
    allocate(this%order,source=order)
    if(size(mat)/=len(abbrev_in)) ERROR STOP "Matrix size has the be the same as the length of abbrev_in"
    allocate(this%dat(size(mat)))
    do i=1,size(mat)
        this%dat(i)%nnz     =mat(i)%nnz
        this%dat(i)%dim_mat =mat(i)%dim_mat
        Call move_alloc(mat(i)%col,this%dat(i)%col)
        Call mat(i)%destroy()
    enddo
    this%N_mode=this%dat(1)%dim_mat(1)
end subroutine
end module
