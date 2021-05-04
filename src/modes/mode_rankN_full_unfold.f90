module m_mode_construction_rankN_full_manual
use m_mode_construction
use m_derived_types, only : lattice,number_different_order_parameters
implicit none
private
public F_mode_rankN_full_manual

type, extends(F_mode) :: F_mode_rankN_full_manual  !contains all entries
    contains
    !necessary routines as defined by class
    procedure   :: get_mode   !subroutine which returns the mode 
    procedure   :: get_mode_exc
    procedure   :: mode_reduce_comp
    procedure   :: get_ind_site

    procedure   :: get_mode_single_cont  !
    procedure   :: get_mode_single_disc

    procedure   :: copy
    procedure   :: destroy
    procedure   :: is_same

    !MPI
    procedure   :: bcast
    procedure   :: send
    procedure   :: recv

    !local construction routine
    procedure   :: init_order
end type

contains

subroutine get_ind_site(this,comp,site,ind)
    class(F_mode_rankN_full_manual),intent(in)  :: this
    integer,intent(in)                          :: comp  !mode index
    integer,intent(in)                          :: site    !entry
    integer,intent(inout),allocatable           :: ind(:)

    integer         :: inner_dim_mode, i

    ERROR STOP "IMPLEMENT"
end subroutine

subroutine get_mode_single_disc(this,lat,comp,site,ind,vec)
    class(F_mode_rankN_full_manual),intent(in)   :: this
    type(lattice),intent(in)                    :: lat
    integer,intent(in)                          :: comp  !mode index
    integer,intent(in)                          :: site    !entry
    integer,intent(inout),allocatable           :: ind(:)
    real(8),intent(inout),allocatable           :: vec(:)

    ERROR STOP "IMPLEMENT"
end subroutine

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


subroutine get_mode_exc(this,lat,comp,vec)
    use, intrinsic :: iso_fortran_env, only : error_unit
    class(F_mode_rankN_full_manual),intent(in)  :: this
    type(lattice),intent(in)                    :: lat       !lattice type which knows about all states
    integer,intent(in)                          :: comp       !which mode is kept
    real(8),intent(inout)                       :: vec(:)

    logical                                     :: exclude(this%N_mode)
    integer                                     :: i

    exclude=.false.
    exclude(comp)=.true.
    Call lat%set_order_comb_exc(this%order,vec,exclude)
end subroutine

subroutine mode_reduce_comp(this,lat,vec_in,comp,vec_out)
    class(F_mode_rankN_full_manual),intent(in)  :: this
    real(8),intent(in)                          :: vec_in(:)
    type(lattice),intent(in)                    :: lat       !lattice type which knows about all states
    integer,intent(in)                          :: comp !of which operator the first entry is kept
    real(8),intent(out)                         :: vec_out(lat%dim_modes(this%order(comp))*lat%Ncell)

    Call lat%reduce(vec_in,this%order,comp,vec_out)
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
    deallocate(this%order)
end subroutine

subroutine copy(this,F_out)
    class(F_mode_rankN_full_manual),intent(in)    :: this
    class(F_mode),allocatable,intent(inout) :: F_out

    Call this%copy_base(F_out) 
    select type(F_out)
    type is(F_mode_rankN_full_manual)
        continue
    class default
        ERROR STOP "FAILED TO COPY F_mode_rankN_full_manualer mode to F_out"
    end select
end subroutine

subroutine init_order(this,lat,abbrev_in)
    use m_derived_types, only: op_abbrev_to_int
    class(F_mode_rankN_full_manual),intent(inout) :: this
    type(lattice),intent(in)                :: lat       !lattice type which knows about all states
    character(len=*), intent(in)            :: abbrev_in
    integer     :: order(len(abbrev_in))
    integer     :: size_in
    integer     :: i

    order=op_abbrev_to_int(abbrev_in)
    size_in=lat%Ncell*product(lat%dim_modes(order))
    Call this%init_base(order,size_in)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!            MPI ROUTINES           !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine bcast(this,comm)
    use mpi_basic                
    class(F_mode_rankN_full_manual),intent(inout) ::  this        !this might fail if the server threads non-allocated class(F_mode), TAKE CARE OF THIS IN HAM_BASE
    type(mpi_type),intent(in)               ::  comm
#ifdef CPP_MPI
    integer     :: ierr
  
    Call this%bcast_base(comm)
#else
    continue
#endif
end subroutine 

subroutine send(this,ithread,tag,com)
!    use mpi_basic                
    class(F_mode_rankN_full_manual),intent(in)    :: this
    integer,intent(in)          :: ithread
    integer,intent(in)          :: tag
    integer,intent(in)          :: com

#ifdef CPP_MPI
    integer     :: ierr

    Call this%send_base(ithread,tag,com)
    ERROR STOP "IMPLEMENT"
#else
    continue
#endif
end subroutine

subroutine recv(this,ithread,tag,com)
!    use mpi_basic                
    class(F_mode_rankN_full_manual),intent(inout) :: this
    integer,intent(in)          :: ithread
    integer,intent(in)          :: tag
    integer,intent(in)          :: com

#ifdef CPP_MPI
    integer     :: ierr

    Call this%recv_base(ithread,tag,com)
    ERROR STOP "IMPLEMENT"
#else
    continue
#endif
end subroutine

end module
