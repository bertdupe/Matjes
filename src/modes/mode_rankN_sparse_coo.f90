module m_mode_construction_rankN_sparse_coo
use m_mode_construction
use m_derived_types, only : lattice,number_different_order_parameters
use m_coo_mat
implicit none
private
public F_mode_rankN_sparse_coo

type, extends(F_mode) :: F_mode_rankN_sparse_coo !contains all entries
    type(coo_mat),allocatable       :: mat(:)
    integer,allocatable             :: order(:)
    integer                         :: N_mode=-1
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
    class(F_mode_rankN_sparse_coo),intent(in)  :: this
    type(lattice),intent(in)                    :: lat       !lattice type which knows about all states
    integer,intent(in)                          :: op_exc !of which operator the first entry is kept
    real(8),intent(inout)                       :: vec(:)

    logical         :: exclude(size(this%order))
    integer         :: N
    real(8)         :: tmp_internal(this%N_mode,size(this%mat))
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

    tmp_internal=0.d0
    do i=1,size(this%mat)
        if(exclude(i))then
            tmp_internal(:,i)=1.0d0
        else
            Call lat%set_order_point(this%order(i),mode_base)
            do j=1,this%mat(i)%nnz
                tmp_internal(this%mat(i)%row(j),i)=tmp_internal(this%mat(i)%row(j),i)+this%mat(i)%val(j)*mode_base(this%mat(i)%col(j))
            enddo
        endif
    enddo
    vec=product(tmp_internal,dim=2)
    nullify(mode_base)
end subroutine

subroutine mode_reduce(this,lat,vec_in,op_keep,vec_out)
    use, intrinsic :: iso_fortran_env, only : error_unit
    class(F_mode_rankN_sparse_coo),intent(in)  :: this
    real(8),intent(in)                          :: vec_in(:)
    type(lattice),intent(in)                    :: lat       !lattice type which knows about all states
    integer,intent(in)                          :: op_keep   !of which operator the first entry is kept
    real(8),intent(out)                         :: vec_out(lat%dim_modes(op_keep)*lat%Ncell)

    integer     ::  i_order,i

    i_order=findloc(this%order,op_keep,dim=1)
    vec_out=0.0d0
    do i=1,this%mat(i_order)%nnz
        vec_out(this%mat(i_order)%col(i))=vec_out(this%mat(i_order)%col(i))+vec_in(this%mat(i_order)%row(i))
    enddo
end subroutine


subroutine get_mode(this,lat,mode,tmp)
    !bad implementation not intended for production
    class(F_mode_rankN_sparse_coo),intent(in)   :: this
    type(lattice),intent(in)                    :: lat       !lattice type which knows about all states
    real(8),intent(out),pointer                 :: mode(:)   !pointer to required mode
    real(8),allocatable,target,intent(inout)    :: tmp(:)

    integer         :: N
    real(8)         :: tmp_internal(this%N_mode,size(this%mat))
    real(8),pointer :: mode_base(:)

    integer         :: i,j

    allocate(tmp(this%N_mode),source=0.0d0)
    mode=>tmp

    tmp_internal=0.d0
    do i=1,size(this%mat)
        Call lat%set_order_point(this%order(i),mode_base)
        do j=1,this%mat(i)%nnz
            tmp_internal(this%mat(i)%row(j),i)=tmp_internal(this%mat(i)%row(j),i)+this%mat(i)%val(j)*mode_base(this%mat(i)%col(j))
        enddo
    enddo
    mode=product(tmp_internal,dim=2)
    nullify(mode_base)
end subroutine

function is_same(this,comp)result(same)
    class(F_mode_rankN_sparse_coo),intent(in)       :: this
    class(F_mode),intent(in)                   :: comp
    logical                                    :: same

    ERROR STOP "IMPLEMENT"
    !same=.false.
    !select type(comp) 
    !type is(F_mode_rankN_sparse_coo)
    !    same=all(this%order==comp%order)
    !end select
end function

subroutine destroy(this)
    class(F_mode_rankN_sparse_coo),intent(inout) ::  this
    integer ::  i
    this%N_mode=-1
    do i=1,size(this%mat)
        Call this%mat(i)%destroy()
    enddo
    deallocate(this%mat)
    deallocate(this%order)
end subroutine

subroutine copy(this,F_out)
    class(F_mode_rankN_sparse_coo),intent(in)    :: this
    class(F_mode),allocatable,intent(inout) :: F_out

    integer ::  i

    if(.not.allocated(F_out)) allocate(F_mode_rankN_sparse_coo::F_out)
    select type(F_out)
    class is(F_mode_rankN_sparse_coo)
        if(.not.allocated(F_out%order)) allocate(F_out%order(size(this%order)))
        if(size(F_out%order)/=size(this%order)) ERROR STOP "CANNOT COPY ORDER AS RANKS DIFFER"
        F_out%order=this%order
        F_out%N_mode=this%N_mode
        allocate(F_out%mat(size(this%mat)))
        do i=1,size(this%mat)
            Call this%mat(i)%copy(F_out%mat(i))
        enddo
    class default
        ERROR STOP "FAILED TO COPY F_mode_rankN_sparse_coo mode to F_out"
    end select
end subroutine

subroutine bcast(this,comm)
    use mpi_basic                
    class(F_mode_rankN_sparse_coo),intent(inout) ::  this        !this might fail if the server threads non-allocated class(F_mode), TAKE CARE OF THIS IN HAM_BASE
    type(mpi_type),intent(in)               ::  comm
#ifdef CPP_MPI
    integer     :: ierr
    integer     ::  N
  
    !THIS MIGHT BE INSUFFICIENT, MAYBE ONE HAS TO CHECK IF THE F_MODE IS ALREADY ALLOCATED TO THE F_mode_rankN_sparse_coo type
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
        Call MPI_Bcast(this%mat(i)%nnz    ,1, MPI_INTEGER, comm%mas, comm%com,ierr)
        Call MPI_Bcast(this%mat(i)%dim_mat,2, MPI_INTEGER, comm%mas, comm%com,ierr)
        if(.not.comm%ismas)then
            allocate(this%mat(i)%row(this%mat(i)%nnz))
            allocate(this%mat(i)%col(this%mat(i)%nnz))
            allocate(this%mat(i)%val(this%mat(i)%nnz))
        endif
        Call MPI_Bcast(this%mat(i)%row,this%mat(i)%nzz, MPI_INTEGER         , comm%mas, comm%com,ierr)
        Call MPI_Bcast(this%mat(i)%col,this%mat(i)%nzz, MPI_INTEGER         , comm%mas, comm%com,ierr)
        Call MPI_Bcast(this%mat(i)%val,this%mat(i)%nzz, MPI_DOUBLE_PRECISION, comm%mas, comm%com,ierr)
    enddo
#else
    continue
#endif
end subroutine 

subroutine init_order(this,lat,abbrev_in,mat)
    use m_derived_types, only: op_abbrev_to_int
    class(F_mode_rankN_sparse_coo),intent(inout) :: this
    type(lattice),intent(in)                :: lat       !lattice type which knows about all states
    character(len=*), intent(in)            :: abbrev_in !considered order abbreviations
    type(coo_mat),intent(inout)             :: mat(:)    !input matrices, destroyed when returned
    integer     :: order(len(abbrev_in))
    integer     :: i     

    order=op_abbrev_to_int(abbrev_in)
    allocate(this%order,source=order)
    if(size(mat)/=len(abbrev_in)) ERROR STOP "Matrix size has the be the same as the length of abbrev_in"
    allocate(this%mat(size(mat)))
    do i=1,size(mat)
        Call mat(i)%mv(this%mat(i))
    enddo
    this%N_mode=this%mat(1)%dim_mat(1)  !OR 2, think about it
end subroutine
end module
