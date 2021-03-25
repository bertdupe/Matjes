module m_mode_construction
! module that contains functions which combines the different sites on the left and the right of the Hamiltonians
! to transform Hamiltonians of ranks >2 to Hamiltonians of rank 2
use m_derived_types, only : lattice, number_different_order_parameters
implicit none

private
public F_mode

type, abstract :: F_mode
    integer             :: order_occ(number_different_order_parameters)=0
    integer,allocatable :: order(:) !order parameter indices of respective rank
    integer             :: N_mode   !size(order) number of order parameters that make up the mode
    contains
    procedure(int_get_mode),deferred            :: get_mode           !subroutine which returns the mode 
    procedure(int_get_mode_exc_ind),deferred    :: get_mode_exc_ind   !subroutine which returns the mode excluding one order parameter 
    procedure(int_mode_reduce_ind),deferred     :: mode_reduce_ind    !subroutine which reduces an input vector to the shape of a single constituent state according the the F_mode rules

    procedure(int_get_mode_single_cont),deferred :: get_mode_single_cont       !subroutine which returns the mode 

    procedure(int_is_same),deferred         :: is_same
    procedure(int_destroy),deferred         :: destroy
    procedure(int_copy),deferred            :: copy
    procedure(int_bcast),deferred           :: bcast

    procedure,NON_OVERRIDABLE               :: get_mode_exc   !subroutine which returns the mode excluding one order parameter 
    procedure,NON_OVERRIDABLE               :: mode_reduce    !subroutine which reduces an input vector to the shape of a single constituent state according the the F_mode rules

    procedure                               :: init_base
    procedure                               :: copy_base

    procedure                               :: reduce_other_exc
end type

abstract interface
    subroutine int_get_mode_single_cont(this,lat,order,i,modes,vec,bnd)
        import F_mode,lattice
        class(F_mode),intent(in)                    :: this
        type(lattice),intent(in)                    :: lat
        integer,intent(in)                          :: order
        integer,intent(in)                          :: i
        real(8),pointer,intent(out)                 :: modes(:)
        integer,intent(out)                         :: bnd(2)
        real(8),allocatable,target,intent(out)      :: vec(:)   !space to allocate array if not single operator
    end subroutine

    subroutine int_get_mode(this,lat,mode,tmp)
        !subroutine which return a pointer to the mode, either saved in the tmp array or at some other independent array (depending on the implementation)
        import F_mode,lattice
        class(F_mode),intent(in)                   :: this
        type(lattice),intent(in)                   :: lat       !lattice type which knows about all states
        real(8),intent(out),pointer                :: mode(:)   !pointer to required mode
        real(8),allocatable,target,intent(inout)   :: tmp(:)    !possible temporary storage for mode pointer
    end subroutine

    subroutine int_get_mode_exc_ind(this,lat,ind,vec)
        !subroutine which gets the mode, excluding the first occation of the state indexed by op_exc
        import F_mode,lattice
        class(F_mode),intent(in)                    :: this
        type(lattice),intent(in)                    :: lat    !lattice type which knows about all states
        integer,intent(in)                          :: ind    !operator index which is not multiplied [1,this%N_mode]
        real(8),intent(inout)                       :: vec(:) !result mode excluding the first state with op_exc
    end subroutine

    subroutine int_mode_reduce_ind(this,lat,vec_in,ind,vec_out)
        !subroutine which reduces the input mode according to the rules of this F_modes construction rules
        ! to be in the basis of order index by ind
        import F_mode,lattice
        class(F_mode),intent(in)        :: this
        real(8),intent(in)              :: vec_in(:)
        type(lattice),intent(in)        :: lat       !lattice type which knows about all states
        integer,intent(in)              :: ind
        real(8),intent(out)             :: vec_out(lat%dim_modes(this%order(ind))*lat%Ncell)
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
end interface

contains


subroutine bcast_base(this,comm)
    use mpi_basic                
    class(F_mode),intent(inout) ::  this        !this might fail if the server threads non-allocated class(F_mode), TAKE CARE OF THIS IN HAM_BASE
    type(mpi_type),intent(in)   ::  comm
#ifdef CPP_MPI
    integer     :: ierr
    integer     :: N
  
    !THIS MIGHT BE INSUFFICIENT, MAYBE ONE HAS TO CHECK IF THE F_MODE IS ALREADY ALLOCATED TO THE F_mode_rankN_full_manual type
    STOP "CHECK IF THIS WORKS WITHOUT PREVIOUS ALLOCATION./type stuff/, on non-master threads"
    if(comm%ismas)then
        if(.not.allocated(this%order)) ERROR STOP "CANNOT BCAST SINCE MASTER ORDER IS NOT ALLOCATED"
    endif
    Call MPI_Bcast(this%N_mode,1, MPI_INTEGER, comm%mas, comm%com,ierr)
    if(ierr/=0) ERROR STOP "MPI BCAST FAILED"
    if(.not.allocated(this%order)) allocate(this%order(this%N_mode))
    Call MPI_Bcast(this%order,N, MPI_INTEGER, comm%mas, comm%com,ierr)
    if(ierr/=0) ERROR STOP "MPI BCAST FAILED"
    Call MPI_Bcast(this%order_occ,number_different_order_parameters, MPI_INTEGER, comm%mas, comm%com,ierr)
    if(ierr/=0) ERROR STOP "MPI BCAST FAILED"
#else
    continue
#endif
end subroutine 

subroutine init_base(this,order)
    class(F_mode),intent(inout) :: this
    integer,intent(in)          :: order(:)
    
    integer ::  i
    allocate(this%order,source=order)
    do i=1,number_different_order_parameters
        this%order_occ(i)=count(this%order==i)
    enddo
    this%N_mode=size(this%order)
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
    integer         :: ind

    ind=findloc(this%order,op_keep,dim=1)
#ifdef CPP_DEBUG
    if(ind<1.or.ind>size(this%order))then
        write(error_unit,'(//A,I6)') "Tried to get mode excluding order no.:", op_exc
        write(error_unit,*) "But the mode only contains the order:", this%order
        ERROR STOP "This makes no sense and should probably prevented earlier in the code"
    endif
#endif
    Call this%get_mode_exc_ind(lat,ind,tmp)
    tmp=vec_other*tmp
    Call this%mode_reduce_ind(lat,tmp,ind,res)
    res=res*real(this%order_occ(op_keep),8)
end subroutine

subroutine get_mode_exc(this,lat,op_exc,vec)
    use, intrinsic :: iso_fortran_env, only : error_unit
    class(F_mode),intent(in)        :: this
    type(lattice),intent(in)        :: lat       !lattice type which knows about all states
    integer,intent(in)              :: op_exc !of which operator the first entry is kept
    real(8),intent(inout)           :: vec(:)

    integer         :: ind

    ind=findloc(this%order,op_exc,dim=1)
#ifdef CPP_DEBUG
    if(size(vec)/=this%mode_size) STOP "mode exc call has wrong size for vector"
    if(ind<1.or.ind>size(this%order))then
        write(error_unit,'(//A,I6)') "Tried to get mode excluding order no.:", op_exc
        write(error_unit,*) "But the mode only contains the order:", this%order
        ERROR STOP "This makes no sense and should probably prevented earlier in the code"
    endif
#endif
    Call this%get_mode_exc_ind(lat,ind,vec)
end subroutine

subroutine mode_reduce(this,lat,vec_in,op_keep,vec_out)
    !reduce mode by first occurance of operato op_keep in mode
    use, intrinsic :: iso_fortran_env, only : error_unit
    class(F_mode),intent(in)        :: this
    real(8),intent(in)              :: vec_in(:)
    type(lattice),intent(in)        :: lat       !lattice type which knows about all states
    integer,intent(in)              :: op_keep   !of which operator the first entry is kept
    real(8),intent(out)             :: vec_out(lat%dim_modes(op_keep)*lat%Ncell)

    integer     ::  ind,i

    ind=findloc(this%order,op_keep,dim=1)
#ifdef CPP_DEBUG
    if(ind<1.or.ind>size(this%order))then
        write(error_unit,'(//A,I6)') "Tried to reduce mode excluding order no.:", op_exc
        write(error_unit,*) "But the mode only contains the order:", this%order
        ERROR STOP "This makes no sense and should probably prevented earlier in the code"
    endif
#endif
    Call this%mode_reduce_ind(lat,vec_in,ind,vec_out)
end subroutine
end module 
