module m_mode_construction_rankN_full_manual
use m_mode_construction
use m_derived_types, only : lattice,number_different_order_parameters
implicit none
private
public F_mode_rankN_full_manual

type, extends(F_mode) :: F_mode_rankN_full_manual  !contains all entries
    integer,allocatable,private     :: order(:)
    integer                         :: N_mode=-1
    contains
    !necessary routines as defined by class
    procedure   :: get_mode   !subroutine which returns the mode 
    procedure   :: get_mode_exc
    procedure   :: mode_reduce  

    procedure   :: get_mode_single_cont  !

    procedure   :: copy
    procedure   :: bcast
    procedure   :: destroy
    procedure   :: is_same
    !local construction routine
    procedure   :: init_order
end type

contains

subroutine get_mode_single_cont(this,lat,order,i,modes,vec,bnd)
    class(F_mode_rankN_full_manual),intent(in)  :: this
    type(lattice),intent(in)                    :: lat
    integer,intent(in)                          :: order
    integer,intent(in)                          :: i
    real(8),pointer,intent(out)                 :: modes(:)
    integer,intent(out)                         :: bnd(2)
    real(8),allocatable,target,intent(out)      :: vec(:)   !space to allocate array if not single operator

    ERROR STOP "IMPLEMENT"
end subroutine


subroutine get_mode_exc(this,lat,op_exc,vec)
    use, intrinsic :: iso_fortran_env, only : error_unit
    class(F_mode_rankN_full_manual),intent(in)  :: this
    type(lattice),intent(in)                    :: lat       !lattice type which knows about all states
    integer,intent(in)                          :: op_exc !of which operator the first entry is kept
    real(8),intent(inout)                       :: vec(:)

    logical                                     :: exclude(size(this%order))
    integer                                     :: i

    if(size(vec)/=this%N_mode) STOP "mode exc call has wrong size for vector"
    exclude=.false.
    i=findloc(this%order,op_exc,dim=1)
    if(i<1.or.i>size(this%order))then
        write(error_unit,'(//A,I6)') "Tried to get mode excluding order no.:", op_exc
        write(error_unit,*) "But the mode only contains the order:", this%order
        ERROR STOP "This makes no sense and should probably prevented earlier in the code"
    endif
    exclude(i)=.true.
    Call lat%set_order_comb_exc(this%order,vec,exclude)
end subroutine

subroutine mode_reduce(this,lat,vec_in,op_keep,vec_out)
    use, intrinsic :: iso_fortran_env, only : error_unit
    class(F_mode_rankN_full_manual),intent(in)  :: this
    real(8),intent(in)                          :: vec_in(:)
    type(lattice),intent(in)                    :: lat       !lattice type which knows about all states
    integer,intent(in)                          :: op_keep   !of which operator the first entry is kept
    real(8),intent(out)                         :: vec_out(lat%dim_modes(op_keep)*lat%Ncell)

    Call lat%reduce(vec_in,this%order,op_keep,vec_out)
end subroutine


subroutine get_mode(this,lat,mode,tmp)
    class(F_mode_rankN_full_manual),intent(in)  :: this
    type(lattice),intent(in)                    :: lat       !lattice type which knows about all states
    real(8),intent(out),pointer                 :: mode(:)   !pointer to required mode
    real(8),allocatable,target,intent(inout)    :: tmp(:)

    integer         ::  N
    allocate(tmp(this%N_mode))
    mode=>tmp
    Call lat%set_order_comb(this%order,mode)
end subroutine

function is_same(this,comp)result(same)
    class(F_mode_rankN_full_manual),intent(in)       :: this
    class(F_mode),intent(in)                   :: comp
    logical                                    :: same

    same=.false.
    select type(comp) 
    type is(F_mode_rankN_full_manual)
        same=all(this%order==comp%order)
    end select
end function

subroutine destroy(this)
    class(F_mode_rankN_full_manual),intent(inout) ::  this
    this%N_mode=-1
    deallocate(this%order)
end subroutine

subroutine copy(this,F_out)
    class(F_mode_rankN_full_manual),intent(in)    :: this
    class(F_mode),allocatable,intent(inout) :: F_out

    Call this%copy_base(F_out) 
    select type(F_out)
    type is(F_mode_rankN_full_manual)
        if(.not.allocated(F_out%order)) allocate(F_out%order(size(this%order)))
        if(size(F_out%order)/=size(this%order)) ERROR STOP "CANNOT COPY ORDER AS RANKS DIFFER"
        F_out%order=this%order
        F_out%N_mode=this%N_mode
    class default
        ERROR STOP "FAILED TO COPY F_mode_rankN_full_manualer mode to F_out"
    end select
end subroutine

subroutine bcast(this,comm)
    use mpi_basic                
    class(F_mode_rankN_full_manual),intent(inout) ::  this        !this might fail if the server threads non-allocated class(F_mode), TAKE CARE OF THIS IN HAM_BASE
    type(mpi_type),intent(in)               ::  comm
#ifdef CPP_MPI
    integer     :: ierr
    integer     ::  N
  
    !THIS MIGHT BE INSUFFICIENT, MAYBE ONE HAS TO CHECK IF THE F_MODE IS ALREADY ALLOCATED TO THE F_mode_rankN_full_manual type
    STOP "CHECK IF THIS WORKS WITHOUT PREVIOUS ALLOCATION./type stuff/, on non-master threads"
    if(comm%ismas)then
        if(.not.allocated(this%order)) ERROR STOP "CANNOT BCAST SINCE MASTER ORDER IS NOT ALLOCATED"
        N=size(this%order)
    endif
    Call MPI_Bcast(N,1, MPI_INTEGER, comm%mas, comm%com,ierr)
    if(.not.allocated(this%order)) allocate(this%order(N))
    Call MPI_Bcast(this%order,N, MPI_INTEGER, comm%mas, comm%com,ierr)
    Call MPI_Bcast(this%N_mode,1, MPI_INTEGER, comm%mas, comm%com,ierr)
    if(ierr/=0) ERROR STOP "MPI BCAST FAILED"
#else
    continue
#endif
end subroutine 

subroutine init_order(this,lat,abbrev_in)
    use m_derived_types, only: op_abbrev_to_int
    class(F_mode_rankN_full_manual),intent(inout) :: this
    type(lattice),intent(in)                :: lat       !lattice type which knows about all states
    character(len=*), intent(in)            :: abbrev_in
    integer     :: order(len(abbrev_in))
    integer     :: i

    integer     :: order_occ(number_different_order_parameters)

    order=op_abbrev_to_int(abbrev_in)
    allocate(this%order,source=order)
    do i=1,number_different_order_parameters
        order_occ(i)=count(order==i)
    enddo
    Call this%init_base(order_occ)
    this%N_mode=lat%Ncell*product(lat%dim_modes(this%order))
end subroutine
end module
