module m_mode_construction_rankN_full_manual
use m_mode_construction
use m_derived_types, only : lattice,number_different_order_parameters
implicit none
private
public F_mode_rankN_full_manual

type, extends(F_mode) :: F_mode_rankN_full_manual  !contains all entries
    integer                         :: mode_size=-1
    contains
    !necessary routines as defined by class
    procedure   :: get_mode   !subroutine which returns the mode 
    procedure   :: get_mode_exc_ind
    procedure   :: mode_reduce_ind

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


subroutine get_mode_exc_ind(this,lat,ind,vec)
    use, intrinsic :: iso_fortran_env, only : error_unit
    class(F_mode_rankN_full_manual),intent(in)  :: this
    type(lattice),intent(in)                    :: lat       !lattice type which knows about all states
    integer,intent(in)                          :: ind       !which mode is kept
    real(8),intent(inout)                       :: vec(:)

    logical                                     :: exclude(this%N_mode)
    integer                                     :: i

    exclude=.false.
    exclude(ind)=.true.
    Call lat%set_order_comb_exc(this%order,vec,exclude)
end subroutine

subroutine mode_reduce_ind(this,lat,vec_in,ind,vec_out)
    class(F_mode_rankN_full_manual),intent(in)  :: this
    real(8),intent(in)                          :: vec_in(:)
    type(lattice),intent(in)                    :: lat       !lattice type which knows about all states
    integer,intent(in)                          :: ind !of which operator the first entry is kept
    real(8),intent(out)                         :: vec_out(lat%dim_modes(this%order(ind))*lat%Ncell)

    Call lat%reduce(vec_in,this%order,ind,vec_out)
end subroutine


subroutine get_mode(this,lat,mode,tmp)
    class(F_mode_rankN_full_manual),intent(in)  :: this
    type(lattice),intent(in)                    :: lat       !lattice type which knows about all states
    real(8),intent(out),pointer                 :: mode(:)   !pointer to required mode
    real(8),allocatable,target,intent(inout)    :: tmp(:)

    integer         ::  N
    allocate(tmp(this%mode_size))
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
    this%mode_size=-1
    deallocate(this%order)
end subroutine

subroutine copy(this,F_out)
    class(F_mode_rankN_full_manual),intent(in)    :: this
    class(F_mode),allocatable,intent(inout) :: F_out

    Call this%copy_base(F_out) 
    select type(F_out)
    type is(F_mode_rankN_full_manual)
        F_out%mode_size=this%mode_size
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
  
    Call bcast_base(this,comm)
    Call MPI_Bcast(this%mode_size,1, MPI_INTEGER, comm%mas, comm%com,ierr)
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

    order=op_abbrev_to_int(abbrev_in)
    Call this%init_base(order)
    this%mode_size=lat%Ncell*product(lat%dim_modes(this%order))
end subroutine
end module
