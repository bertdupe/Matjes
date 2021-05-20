module m_mode_construction
! module that contains functions which combines the different sites on the left and the right of the Hamiltonians
! to transform Hamiltonians of ranks >2 to Hamiltonians of rank 2
use m_derived_types, only : lattice, number_different_order_parameters
use, intrinsic :: iso_fortran_env, only : error_unit
implicit none

private
public F_mode

type, abstract :: F_mode
    integer             :: order_occ(number_different_order_parameters)=0
    integer             :: N_mode=-1   !size(order) number of order parameters that make up the mode
    integer             :: mode_size=-1 !overall size of this mode
    integer,allocatable :: order(:) !order parameter indices of respective rank
    contains
    procedure(int_get_mode),deferred            :: get_mode             !subroutine which returns the mode 
    procedure(int_get_mode_exc),deferred        :: get_mode_exc         !subroutine which returns the mode excluding one order parameter 
    procedure(int_mode_reduce_comp),deferred    :: mode_reduce_comp     !subroutine which reduces an input vector to the shape of a single constituent state according the the F_mode rules
    procedure(int_get_ind_site),deferred        :: get_ind_site 
    procedure(int_get_ind_site_expl),deferred   :: get_ind_site_expl

    procedure(int_get_mode_single),deferred      :: get_mode_single  !subroutine which returns the mode as discontiguous array
    procedure,NON_OVERRIDABLE                    :: get_mode_single_expl  !returns discontiguous array of mode to explicit output arrays
    procedure(int_get_mode_single_size),deferred :: get_mode_single_size  !get size necessary to set a single mode

!subroutines resolving order to first index
    procedure,NON_OVERRIDABLE               :: get_mode_exc_op   !subroutine which returns the mode excluding one order parameter 
    procedure,NON_OVERRIDABLE               :: mode_reduce      !subroutine which reduces an input vector to the shape of a single constituent state according the the F_mode rules

    procedure                               :: reduce_other_exc
!routines, which should be implemented more efficiently locally
    procedure                               :: get_mode_disc
    procedure                               :: get_mode_disc_expl
    procedure                               :: get_mode_exc_op_disc
    procedure                               :: get_mode_exc_disc
    procedure                               :: mode_reduce_comp_disc
!misc
    procedure                               :: ind_site
    procedure                               :: get_comp
!basic functionalities
    procedure(int_is_same),deferred         :: is_same
    procedure(int_destroy),deferred         :: destroy
    procedure(int_copy),deferred            :: copy
    procedure(int_bcast),deferred           :: bcast
    procedure                               :: init_base
    procedure                               :: copy_base


!MPI
    procedure(int_send),deferred            :: send
    procedure(int_recv),deferred            :: recv
    procedure                               :: bcast_base
    procedure                               :: send_base
    procedure                               :: recv_base
end type

abstract interface

    subroutine int_get_ind_site(this,comp,site,ind)
        import F_mode
        class(F_mode),intent(in)                    :: this
        integer,intent(in)                          :: comp  !mode component
        integer,intent(in)                          :: site    !entry
        integer,intent(inout),allocatable           :: ind(:)
    end subroutine

    subroutine int_get_ind_site_expl(this,comp,site,size_out,ind)
        import F_mode
        class(F_mode),intent(in)                    :: this
        integer,intent(in)                          :: comp  !mode component
        integer,intent(in)                          :: site    !entry
        integer,intent(in)                          :: size_out
        integer,intent(out)                         :: ind(size_out)
    end subroutine


    subroutine int_get_mode_single(this,lat,comp,site,ind,vec)
        import F_mode, lattice
        class(F_mode),intent(in)                :: this
        type(lattice),intent(in)                :: lat
        integer,intent(in)                      :: comp  !mode index
        integer,intent(in)                      :: site    !entry
        integer,intent(inout),allocatable       :: ind(:)
        real(8),intent(inout),allocatable       :: vec(:)
    end subroutine

    subroutine int_get_mode_single_size(this,order,dim_mode)
        import F_mode
        class(F_mode),intent(in)    :: this
        integer,intent(in)          :: order
        integer,intent(out)         :: dim_mode
    end subroutine

    subroutine int_get_mode(this,lat,mode,tmp)
        !subroutine which return a pointer to the mode, either saved in the tmp array or at some other independent array (depending on the implementation)
        import F_mode,lattice
        class(F_mode),intent(in)                   :: this
        type(lattice),intent(in)                   :: lat       !lattice type which knows about all states
        real(8),intent(out),pointer                :: mode(:)   !pointer to required mode
        real(8),allocatable,target,intent(inout)   :: tmp(:)    !possible temporary storage for mode pointer
    end subroutine

    subroutine int_get_mode_exc(this,lat,comp,vec)
        !subroutine which gets the mode, excluding the first occation of the state indexed by op_exc
        import F_mode,lattice
        class(F_mode),intent(in)                    :: this
        type(lattice),intent(in)                    :: lat    !lattice type which knows about all states
        integer,intent(in)                          :: comp    !operator index which is not multiplied [1,this%N_mode]
        real(8),intent(inout)                       :: vec(:) !result mode excluding the first state with op_exc
    end subroutine

    subroutine int_mode_reduce_comp(this,lat,vec_in,comp,vec_out)
        !subroutine which reduces the input mode according to the rules of this F_modes construction rules
        ! to be in the basis of order index by comp
        import F_mode,lattice
        class(F_mode),intent(in)        :: this
        real(8),intent(in)              :: vec_in(:)
        type(lattice),intent(in)        :: lat       !lattice type which knows about all states
        integer,intent(in)              :: comp
        real(8),intent(out)             :: vec_out(lat%dim_modes(this%order(comp))*lat%Ncell)
    end subroutine

    function int_is_same(this,comp)result(same)
        import F_mode
        class(F_mode),intent(in)                   :: this
        class(F_mode),intent(in)                   :: comp
        logical                                    :: same
    end function 

    subroutine int_bcast(this,comm)
        use mpi_basic                
        import F_mode
        class(F_mode),intent(inout) ::  this
        type(mpi_type),intent(in)   ::  comm
    end subroutine
    
    subroutine int_copy(this,F_out)
        import F_mode
        class(F_mode),intent(in)                :: this
        class(F_mode),intent(inout),allocatable :: F_out
    end subroutine

    subroutine int_destroy(this)
        import F_mode
        class(F_mode),intent(inout) ::  this
    end subroutine

    subroutine int_send(this,ithread,tag,com)
        import F_mode
        class(F_mode),intent(in)    :: this
        integer,intent(in)          :: ithread
        integer,intent(in)          :: tag
        integer,intent(in)          :: com
    end subroutine

    subroutine int_recv(this,ithread,tag,com)
        import F_mode
        class(F_mode),intent(inout) :: this
        integer,intent(in)          :: ithread
        integer,intent(in)          :: tag
        integer,intent(in)          :: com
    end subroutine

end interface

contains


subroutine bcast_base(this,comm)
    use mpi_basic                
    class(F_mode),intent(inout) ::  this        !this might fail if the server threads non-allocated class(F_mode), TAKE CARE OF THIS IN HAM_BASE
    type(mpi_type),intent(in)   ::  comm
#ifdef CPP_MPI
    integer     :: ierr
 
    if(comm%ismas)then
        if(.not.allocated(this%order)) ERROR STOP "CANNOT BCAST SINCE MASTER ORDER IS NOT ALLOCATED"
    endif
    Call MPI_Bcast(this%order_occ, number_different_order_parameters, MPI_INTEGER, comm%mas, comm%com,ierr)
    if(ierr/=0) ERROR STOP "MPI BCAST FAILED"
    Call MPI_Bcast(this%mode_size, 1, MPI_INTEGER, comm%mas, comm%com,ierr)
    if(ierr/=0) ERROR STOP "MPI BCAST FAILED"
    Call MPI_Bcast(this%N_mode   , 1, MPI_INTEGER, comm%mas, comm%com,ierr)
    if(ierr/=0) ERROR STOP "MPI BCAST FAILED"

    if(.not.allocated(this%order)) allocate(this%order(this%N_mode))
    Call MPI_Bcast(this%order, this%N_mode, MPI_INTEGER, comm%mas, comm%com,ierr)
    if(ierr/=0) ERROR STOP "MPI BCAST FAILED"
#else
    continue
#endif
end subroutine 

subroutine init_base(this,order,size_in)
    class(F_mode),intent(inout) :: this
    integer,intent(in)          :: order(:)
    integer,intent(in)          :: size_in
    
    integer ::  i
    allocate(this%order,source=order)
    do i=1,number_different_order_parameters
        this%order_occ(i)=count(this%order==i)
    enddo
    this%N_mode=size(this%order)
    this%mode_size=size_in
end subroutine

subroutine copy_base(this,F_out)
    class(F_mode),intent(in)                :: this
    class(F_mode),intent(inout),allocatable :: F_out
    
    if(.not.allocated(F_out)) allocate(F_out,mold=this)
    if(.not.allocated(F_out%order)) allocate(F_out%order(size(this%order)))
    if(size(F_out%order)/=size(this%order)) ERROR STOP "CANNOT COPY ORDER AS RANKS DIFFER"
    F_out%order=this%order
    F_out%order_occ=this%order_occ
    F_out%N_mode=this%N_mode
    F_out%mode_size=this%mode_size
end subroutine

subroutine reduce_other_exc(this,lat,op_keep,vec_other,res)
    !subroutine to get the derivative, may be overwritten
    class(F_mode),intent(in)        :: this
    type(lattice),intent(in)        :: lat          !lattice type which knows about all states
    integer,intent(in)              :: op_keep      !operator index which is to be kept (1,number_different_order_parameters)
    real(8),intent(in)              :: vec_other(:) !other vector to which the Hamiltonian has been multiplied
    real(8),intent(inout)           :: res(:)

    real(8)         :: tmp(size(vec_other))
    integer         :: i
    integer         :: comp

    comp=findloc(this%order,op_keep,dim=1)
#ifdef CPP_DEBUG
    if(comp<1.or.comp>size(this%order))then
        write(error_unit,'(//A,I6)') "Tried to get mode keeping order no.:", op_keep
        write(error_unit,*) "But the mode only contains the order:", this%order
        ERROR STOP "This makes no sense and should probably prevented earlier in the code"
    endif
#endif
    Call this%get_mode_exc(lat,comp,tmp)
    tmp=vec_other*tmp
    Call this%mode_reduce_comp(lat,tmp,comp,res)
    res=res*real(this%order_occ(op_keep),8)
end subroutine

subroutine get_mode_exc_op(this,lat,op_exc,vec)
    class(F_mode),intent(in)        :: this     
    type(lattice),intent(in)        :: lat      !lattice type which knows about all states
    integer,intent(in)              :: op_exc   !of which operator the first entry is kept
    real(8),intent(inout)           :: vec(:)

    integer         :: comp

    comp=findloc(this%order,op_exc,dim=1)
#ifdef CPP_DEBUG
!    if(size(vec)/=this%mode_size) STOP "mode exc call has wrong size for vector"
    if(comp<1.or.comp>size(this%order))then
        write(error_unit,'(//A,I6)') "Tried to get mode excluding order no.:", op_exc
        write(error_unit,*) "But the mode only contains the order:", this%order
        ERROR STOP "This makes no sense and should probably prevented earlier in the code"
    endif
#endif
    Call this%get_mode_exc(lat,comp,vec)
end subroutine

subroutine get_mode_exc_op_disc(this,lat,op_exc,ind,vec)
    class(F_mode),intent(in)        :: this     
    type(lattice),intent(in)        :: lat      !lattice type which knows about all states
    integer,intent(in)              :: op_exc   !of which operator the first entry is kept
    integer,intent(in)              :: ind(:)
    real(8),intent(inout)           :: vec(:)

    real(8)                         :: tmp(this%mode_size)
    integer         :: comp

    comp=findloc(this%order,op_exc,dim=1)
#ifdef CPP_DEBUG
!    if(size(vec)/=this%mode_size) STOP "mode exc call has wrong size for vector"
    if(comp<1.or.comp>size(this%order))then
        write(error_unit,'(//A,I6)') "Tried to get mode excluding order no.:", op_exc
        write(error_unit,*) "But the mode only contains the order:", this%order
        ERROR STOP "This makes no sense and should probably prevented earlier in the code"
    endif
#endif
    Call this%get_mode_exc(lat,comp,tmp)
    vec=tmp(ind)
end subroutine

subroutine get_mode_exc_disc(this,lat,comp,ind,vec)
    class(F_mode),intent(in)        :: this     
    type(lattice),intent(in)        :: lat      !lattice type which knows about all states
    integer,intent(in)              :: comp
    integer,intent(in)              :: ind(:)
    real(8),intent(inout)           :: vec(:)

    real(8)                         :: tmp(this%mode_size)

    Call this%get_mode_exc(lat,comp,tmp)
    vec=tmp(ind)
end subroutine

subroutine mode_reduce(this,lat,vec_in,op_keep,vec_out)
    !reduce mode by first occurance of operato op_keep in mode
    class(F_mode),intent(in)        :: this
    real(8),intent(in)              :: vec_in(:)
    type(lattice),intent(in)        :: lat       !lattice type which knows about all states
    integer,intent(in)              :: op_keep   !of which operator the first entry is kept
    real(8),intent(out)             :: vec_out(lat%dim_modes(op_keep)*lat%Ncell)

    integer     ::  comp,i

    comp=findloc(this%order,op_keep,dim=1)
#ifdef CPP_DEBUG
    if(comp<1.or.comp>size(this%order))then
        write(error_unit,'(//A,I6)') "Tried to reduce mode keeping order no.:", op_keep
        write(error_unit,*) "But the mode only contains the order:", this%order
        ERROR STOP "This makes no sense and should probably prevented earlier in the code"
    endif
#endif
    Call this%mode_reduce_comp(lat,vec_in,comp,vec_out)
end subroutine

subroutine mode_reduce_comp_disc(this,ind_in,vec_in,comp,ind_out,vec_out)
    class(F_mode),intent(in)    :: this
    integer,intent(in)          :: ind_in(:)
    real(8),intent(in)          :: vec_in(:)
    integer,intent(in)          :: comp 
    integer,intent(in)          :: ind_out(:)
    real(8),intent(inout)       :: vec_out(:)

    ERROR STOP "IMPLEMENT LOCALLY"

end subroutine

subroutine get_mode_disc(this,lat,ind,vec)
    !should be implemented more efficiently for relevant derived types
    class(F_mode),intent(in)                :: this
    type(lattice),intent(in)                :: lat
    integer,intent(in)                      :: ind(:)
    real(8),intent(inout),allocatable       :: vec(:)

    real(8),pointer                         :: mode(:)   !pointer to required mode
    real(8),allocatable,target              :: tmp(:)    !possible temporary storage for mode pointer

    if(.not.allocated(vec)) allocate(vec(size(ind)))
    Call this%get_mode(lat,mode,tmp)
    vec=mode(ind)
end subroutine

subroutine get_mode_disc_expl(this,lat,N,ind,vec)
    !should be implemented more efficiently for relevant derived types
    class(F_mode),intent(in)                :: this
    type(lattice),intent(in)                :: lat
    integer,intent(in)                      :: N
    integer,intent(in)                      :: ind(N)
    real(8),intent(out)                     :: vec(N)

    real(8),pointer                         :: mode(:)   !pointer to required mode
    real(8),allocatable,target              :: tmp(:)    !possible temporary storage for mode pointer

    Call this%get_mode(lat,mode,tmp)
    vec=mode(ind)
end subroutine


subroutine ind_site(this,op_keep,site,ind)
    class(F_mode),intent(in)            :: this
    integer,intent(in)                  :: op_keep  !mode index
    integer,intent(in)                  :: site    !entry
    integer,intent(inout),allocatable   :: ind(:)
    integer     :: comp

    comp=findloc(this%order,op_keep,dim=1)
#ifdef CPP_DEBUG
    if(comp<1.or.comp>size(this%order))then
        write(error_unit,'(//A,I6)') "Tried to reduce mode keeping order no.:", op_keep
        write(error_unit,*) "But the mode only contains the order:", this%order
        ERROR STOP "This makes no sense and should probably prevented earlier in the code"
    endif
#endif
    Call this%get_ind_site(comp,site,ind)
end subroutine

function get_comp(this,op_keep)result(comp)
    class(F_mode),intent(in)            :: this
    integer,intent(in)                  :: op_keep  !mode index
    integer     :: comp

    comp=findloc(this%order,op_keep,dim=1)
#ifdef CPP_DEBUG
    if(comp<1.or.comp>size(this%order))then
        write(error_unit,'(//A,I6)') "Tried to reduce mode keeping order no.:", op_keep
        write(error_unit,*) "But the mode only contains the order:", this%order
        ERROR STOP "This makes no sense and should probably prevented earlier in the code"
    endif
#endif
end function 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!            MPI ROUTINES           !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine send_base(this,ithread,tag,com)
    use mpi_basic                
    class(F_mode),intent(in)    :: this
    integer,intent(in)          :: ithread
    integer,intent(in)          :: tag
    integer,intent(in)          :: com

#ifdef CPP_MPI
    integer     :: ierr
    integer     :: N

    N=size(this%order)
    Call MPI_SEND(this%order_occ,     number_different_order_parameters, MPI_INTEGER,   ithread, tag, com, ierr)
    Call MPI_SEND(this%N_mode,        1, MPI_INTEGER,   ithread, tag, com, ierr)
    Call MPI_SEND(this%mode_size,     1, MPI_INTEGER,   ithread, tag, com, ierr)
    Call MPI_SEND(N,                  1, MPI_INTEGER,   ithread, tag, com, ierr)
    Call MPI_SEND(this%order,         N, MPI_INTEGER,   ithread, tag, com, ierr)
#else
    continue
#endif
end subroutine

subroutine recv_base(this,ithread,tag,com)
    use mpi_basic                
    class(F_mode),intent(inout) :: this
    integer,intent(in)          :: ithread
    integer,intent(in)          :: tag
    integer,intent(in)          :: com

#ifdef CPP_MPI
    integer     :: ierr
    integer     :: N
    integer     :: stat(MPI_STATUS_SIZE) 

    Call MPI_RECV(this%order_occ,     number_different_order_parameters, MPI_INTEGER,   ithread, tag, com, stat, ierr)
    Call MPI_RECV(this%N_mode,        1, MPI_INTEGER,   ithread, tag, com, stat, ierr)
    Call MPI_RECV(this%mode_size,     1, MPI_INTEGER,   ithread, tag, com, stat, ierr)
    Call MPI_RECV(N,                  1, MPI_INTEGER,   ithread, tag, com, stat, ierr)
    allocate(this%order(N))
    Call MPI_RECV(this%order,         N, MPI_INTEGER,   ithread, tag, com, stat, ierr)
#else
    continue
#endif
end subroutine

subroutine get_mode_single_expl(this,lat,comp,site,size_out,ind,vec)
    !returns discontiguous array of mode to explicit output arrays
    class(F_mode),intent(in)    :: this
    type(lattice),intent(in)    :: lat
    integer,intent(in)          :: comp  !mode index
    integer,intent(in)          :: site    !entry
    integer,intent(in)          :: size_out  !inner dim_mode
    integer,intent(out)         :: ind(size_out)
    real(8),intent(out)         :: vec(size_out)

    Call this%get_ind_site_expl(comp,site,size_out,ind)
    Call this%get_mode_disc_expl(lat,size_out,ind,vec)
end subroutine
end module 
