module m_mode_construction_rank1_point
use m_type_lattice, only: dim_modes_inner
use m_mode_construction
use m_derived_types, only : lattice,number_different_order_parameters
use, intrinsic  ::  ISO_FORTRAN_ENV, only: error_unit
implicit none
private
public F_mode_rank1_point

type, extends(F_mode) :: F_mode_rank1_point
    ! uses order(1) from the base type to get the correct order from the lattice type
    contains
    !necessary routines as defined by class
    procedure   :: get_mode   !subroutine which returns the mode 
    procedure   :: mode_reduce_comp
    procedure   :: get_ind_site
    !non-sensical routines required by class (exclude makes no sense to rank1 -> always returns 1)
    procedure   :: get_mode_exc
    procedure   :: get_mode_exc_disc


    procedure   :: get_mode_single_size
    procedure   :: get_mode_disc
    procedure   :: reduce_site_vec

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


subroutine reduce_site_vec(this,comp,vec_in,vec_out)
    !reduce the vector along the comp mode if the vec_in has been set with get_ind_site
    !doesn't really make sense for rank1 since there is nothing to reduce
    class(F_mode_rank1_point),intent(in)        :: this
    integer,intent(in)                          :: comp
    real(8),intent(in)                          :: vec_in(:)
    real(8),intent(out)                         :: vec_out(:)

    vec_out=vec_in      
end subroutine


subroutine get_mode_single_size(this,order,dim_mode)
    !returns the size of the vector necessary to get to mode set by a single site
    class(F_mode_rank1_point),intent(in)    :: this
    integer,intent(in)          :: order
    integer,intent(out)         :: dim_mode

    if(order/=this%order(1))then
#ifdef CPP_DEBUG
        write(error_unit,'(A)') "trying to get single mode of order which is not the order of the F_mode"
#endif
        dim_mode=0
    else
        dim_mode=dim_modes_inner(order)
    endif
end subroutine

subroutine get_mode_disc(this,lat,N,ind,vec)
    use, intrinsic :: iso_c_binding, only: C_PTR, C_loc
    use m_type_lattice, only: dim_modes_inner
    class(F_mode_rank1_point),intent(in)    :: this
    type(lattice),intent(in)                :: lat
    integer,intent(in)                      :: N
    integer,intent(in)                      :: ind(N)
    real(8),intent(out)                     :: vec(N)

    real(8),pointer :: mode_base(:)

    Call lat%set_order_point(this%order(1),mode_base)
    vec=mode_base(ind)
end subroutine

subroutine get_ind_site(this,comp,site,size_out,ind)
    !get the indices corresponding to the 
    class(F_mode_rank1_point),intent(in)    :: this
    integer,intent(in)                      :: comp
    integer,intent(in)                      :: site    !entry
    integer,intent(in)                      :: size_out
    integer,intent(out)                     :: ind(size_out)

    integer         :: i

    ind=(site-1)*size_out+[(i,i=1,size_out)]
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

subroutine get_mode_exc_disc(this,lat,comp,N,ind,vec)
    !doesn't really make sense to call this
    use, intrinsic :: iso_c_binding, only: C_PTR, C_loc
    use m_type_lattice, only: dim_modes_inner
    class(F_mode_rank1_point),intent(in)   :: this
    type(lattice),intent(in)                    :: lat
    integer,intent(in)                          :: comp
    integer,intent(in)                          :: N
    integer,intent(in)                          :: ind(N)
    real(8),intent(out)                         :: vec(N)

#ifdef CPP_DEBUG
    write(error_unit,*) "WARNING, using reduce_site_vec for rank1 mode which is superfluous"
    if(comp/=1) ERROR STOP "COMP shall not be larger than 1 for mode_rank1 get_mode_exc_disc"
#endif
    vec=1.0
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

subroutine init_order(this,lat,abbrev_in)
    use m_derived_types, only: op_abbrev_to_int
    class(F_mode_rank1_point),intent(inout) :: this
    type(lattice),intent(in)                :: lat
    character(len=1), intent(in)            :: abbrev_in
    integer                                 :: order(1)

    order=op_abbrev_to_int(abbrev_in)
    Call this%init_base(order,lat%dim_modes(order(1)))
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!            MPI ROUTINES           !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

subroutine send(this,ithread,tag,com)
!    use mpi_basic                
    class(F_mode_rank1_point),intent(in)    :: this
    integer,intent(in)          :: ithread
    integer,intent(in)          :: tag
    integer,intent(in)          :: com

#ifdef CPP_MPI
    integer     :: ierr

    Call this%send_base(ithread,tag,com)
#else
    continue
#endif
end subroutine

subroutine recv(this,ithread,tag,com)
!    use mpi_basic                
    class(F_mode_rank1_point),intent(inout) :: this
    integer,intent(in)          :: ithread
    integer,intent(in)          :: tag
    integer,intent(in)          :: com

#ifdef CPP_MPI
    integer     :: ierr

    Call this%recv_base(ithread,tag,com)
#else
    continue
#endif
end subroutine



end module
