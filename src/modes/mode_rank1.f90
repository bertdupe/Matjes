module m_mode_construction_rank1_point
use m_type_lattice, only: dim_modes_inner
use m_mode_construction
use m_derived_types, only : lattice,number_different_order_parameters
implicit none
private
public F_mode_rank1_point

type, extends(F_mode) :: F_mode_rank1_point
    contains
    !necessary routines as defined by class
    procedure   :: get_mode   !subroutine which returns the mode 
    procedure   :: get_mode_exc
    procedure   :: mode_reduce_comp
    procedure   :: get_ind_site

    procedure   :: get_mode_single_cont
    procedure   :: get_mode_single_disc

    procedure   :: copy
    procedure   :: bcast
    procedure   :: destroy
    procedure   :: is_same

    !local construction routine
    procedure   :: init_order
end type

contains


subroutine get_ind_site(this,comp,site,ind)
    class(F_mode_rank1_point),intent(in)   :: this
    integer,intent(in)                          :: comp  !mode index
    integer,intent(in)                          :: site    !entry
    integer,intent(inout),allocatable           :: ind(:)

    integer         :: inner_dim_mode, i

    inner_dim_mode=dim_modes_inner(this%order(comp))
    if(.not.allocated(ind)) allocate(ind(inner_dim_mode),source=0)
    ind=[((site-1)*inner_dim_mode+i,i=1,inner_dim_mode)]
end subroutine


subroutine get_mode_single_disc(this,lat,comp,site,ind,vec)
    class(F_mode_rank1_point),intent(in)   :: this
    type(lattice),intent(in)               :: lat
    integer,intent(in)                     :: comp  !mode index
    integer,intent(in)                     :: site    !entry
    integer,intent(inout),allocatable      :: ind(:)
    real(8),intent(inout),allocatable      :: vec(:)

    integer                     :: inner_dim_mode, N_out, i
    real(8),contiguous,pointer  :: mode(:)

#ifdef CPP_DEBUG
    if(comp/=1) ERROR STOP "DOESN'T MAKE SENCE FOR comp /=1"
#endif
    inner_dim_mode=dim_modes_inner(this%order(1))
    if(.not.allocated(ind)) allocate(ind(inner_dim_mode))
    if(.not.allocated(vec)) allocate(vec(inner_dim_mode))
    ind=[((site-1)*inner_dim_mode+i,i=1,inner_dim_mode)]

    Call lat%set_order_point(this%order(1),mode)
    vec=mode(ind)
end subroutine


subroutine get_mode_single_cont(this,lat,order,i,modes,vec,bnd)
    class(F_mode_rank1_point),intent(in)        :: this
    type(lattice),intent(in)                    :: lat
    integer,intent(in)                          :: order
    integer,intent(in)                          :: i
    real(8),pointer,intent(out)                 :: modes(:)
    real(8),allocatable,target,intent(out)      :: vec(:)   !space to allocate array if not single operator
    integer,intent(out)                         :: bnd(2)

#ifdef CPP_DEBUG
    if(order/=this%order(1))then
        ERROR STOP "trying to get single mode of order which is not the order of the F_mode"
    endif
#endif
    Call lat%set_order_point_single_inner(order,i,modes,bnd)
end subroutine

subroutine get_mode(this,lat,mode,tmp)
    class(F_mode_rank1_point),intent(in)       :: this
    type(lattice),intent(in)                   :: lat       !lattice type which knows about all states
    real(8),intent(out),pointer                :: mode(:)   !pointer to required mode
    real(8),allocatable,target,intent(inout)   :: tmp(:)    !not used here

    Call lat%set_order_point(this%order(1),mode)
end subroutine

subroutine get_mode_exc(this,lat,comp,vec)
    use, intrinsic :: iso_fortran_env, only : error_unit
    class(F_mode_rank1_point),intent(in)        :: this
    type(lattice),intent(in)                    :: lat      !lattice type which knows about all states
    integer,intent(in)                          :: comp 
    real(8),intent(inout)                       :: vec(:)

    ERROR STOP "Calling get_mode_ext does not make sense for a rank1 mode"
end subroutine

subroutine mode_reduce_comp(this,lat,vec_in,comp,vec_out)
    use, intrinsic :: iso_fortran_env, only : error_unit
    class(F_mode_rank1_point),intent(in)        :: this
    real(8),intent(in)                          :: vec_in(:)
    type(lattice),intent(in)                    :: lat       !lattice type which knows about all states
    integer,intent(in)                          :: comp      !of which operator the first entry is kept
    real(8),intent(out)                         :: vec_out(lat%dim_modes(this%order(comp))*lat%Ncell)

    ERROR STOP "Calling mode_reduce does not make sense for a rank1 mode"
end subroutine

function is_same(this,comp)result(same)
    class(F_mode_rank1_point),intent(in)       :: this
    class(F_mode),intent(in)                   :: comp
    logical                                    :: same

    same=.false.
    select type(comp) 
    type is(F_mode_rank1_point)
        same=all(this%order==comp%order)
    end select
end function

subroutine destroy(this)
    !nothing really has to be done here
    class(F_mode_rank1_point),intent(inout) ::  this
    this%order=0
end subroutine

subroutine copy(this,F_out)
    class(F_mode_rank1_point),intent(in)    :: this
    class(F_mode),allocatable,intent(inout) :: F_out

    Call this%copy_base(F_out)
    select type(F_out)
    class is(F_mode_rank1_point)
        continue !nothing to do here
    class default
        ERROR STOP "FAILED TO COPY F_mode_rank1_pointer mode to F_out"
    end select
end subroutine

subroutine bcast(this,comm)
    use mpi_basic                
    class(F_mode_rank1_point),intent(inout) ::  this        !this might fail if the server threads non-allocated class(F_mode), TAKE CARE OF THIS IN HAM_BASE
    type(mpi_type),intent(in)               ::  comm
#ifdef CPP_MPI
    integer     :: ierr

    Call this%bcast_base(comm)
#else
    continue
#endif
end subroutine 

subroutine init_order(this,lat,abbrev_in)
    use m_derived_types, only: op_abbrev_to_int
    class(F_mode_rank1_point),intent(inout) :: this
    type(lattice),intent(in)                :: lat
    character(len=1), intent(in)            :: abbrev_in
    integer                                 :: order(1)

    order=op_abbrev_to_int(abbrev_in)
    Call this%init_base(order,lat%dim_modes(order(1)))
end subroutine
end module
