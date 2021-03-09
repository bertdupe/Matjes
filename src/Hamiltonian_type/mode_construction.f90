module m_mode_construction
! module that contains functions which combines the different sites on the left and the right of the Hamiltonians
! to transform Hamiltonians of ranks >2 to Hamiltonians of rank 2
use m_derived_types, only : lattice
implicit none

type, abstract :: F_mode
    contains
    procedure(int_get_mode),deferred          :: get_mode   !subroutine which returns the mode 
    !TODO: something which compares the mode construction when checking of Hamiltonian terms can be combined

    procedure(int_is_same),deferred         :: is_same
    procedure(int_destroy),deferred         :: destroy
    procedure(int_copy),deferred            :: copy
    procedure(int_bcast),deferred           :: bcast

end type

abstract interface
    subroutine int_get_mode(this,lat,mode,tmp)
        import F_mode,lattice
        class(F_mode),intent(in)                   :: this
        type(lattice),intent(in)                   :: lat       !lattice type which knows about all states
        real(8),intent(out),pointer                :: mode(:)   !pointer to required mode
        real(8),allocatable,target,intent(inout)   :: tmp(:)    !possible temporary storage for mode pointer
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

end module 


module m_mode_construction_rank1_point
use m_mode_construction
use m_derived_types, only : lattice,number_different_order_parameters
implicit none
private
public F_mode_rank1_point

type, extends(F_mode) :: F_mode_rank1_point
    integer,private     ::  order=0
    contains
    !necessary routines as defined by class
    procedure   :: get_mode   !subroutine which returns the mode 
    procedure   :: copy
    procedure   :: bcast
    procedure   :: destroy
    procedure   :: is_same

    !local construction routine
    procedure   :: init_order

end type

contains

subroutine get_mode(this,lat,mode,tmp)
    class(F_mode_rank1_point),intent(in)       :: this
    type(lattice),intent(in)                   :: lat       !lattice type which knows about all states
    real(8),intent(out),pointer                :: mode(:)   !pointer to required mode
    real(8),allocatable,target,intent(inout)   :: tmp(:)    !not used here

    Call lat%set_order_point(this%order,mode)
end subroutine

function is_same(this,comp)result(same)
    class(F_mode_rank1_point),intent(in)       :: this
    class(F_mode),intent(in)                   :: comp
    logical                                    :: same

    same=.false.
    select type(comp) 
    type is(F_mode_rank1_point)
        same=this%order==comp%order
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

    if(.not.allocated(F_out)) allocate(F_mode_rank1_point::F_out)
    select type(F_out)
    type is(F_mode_rank1_point)
        F_out%order=this%order
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
  
!    if(comm%is_mas)then
!        if(.not.allocated(this)) ERROR STOP "Cannot bcast F_mode_rank1_point if the master entry is not allocated"
!    else
!        if(.not.allocated(this)) allocate(F_mode_rank1_point::this)
!    endif
    !THIS MIGHT BE INSUFFICIENT, MAYBE ONE HAS TO CHECK IF THE F_MODE IS ALREADY ALLOCATED TO THE F_mode_rank1_point type
    STOP "CHECK IF THIS WORKS WITHOUT PREVIOUS ALLOCATION./type stuff/, on non-master threads"
    Call MPI_Bcast(this%order,1, MPI_INTEGER, comm%mas, comm%com,ierr)
    if(ierr/=0) ERROR STOP "MPI BCAST FAILED"
#else
    continue
#endif
end subroutine 

subroutine init_order(this,abbrev_in)
    use m_derived_types, only: op_abbrev_to_int
    class(F_mode_rank1_point),intent(inout) :: this
    character(len=1), intent(in)            :: abbrev_in
    integer                                 :: order(1)

    order=op_abbrev_to_int(abbrev_in)
    this%order=order(1)
end subroutine
end module

module m_mode_public
use m_mode_construction, only: F_mode
use m_mode_construction_rank1_point, only: F_mode_rank1_point
public

contains 
subroutine mode_set_rank1(mode,abbrev_in)
    !function to set mode to single rank avoiding ecasing all the allocation and type select stuff to be more concise in the Energy definitions
    use, intrinsic :: iso_fortran_env, only : error_unit
    class(F_mode),allocatable,intent(inout) :: mode
    character(len=1), intent(in)            :: abbrev_in

    if(allocated(mode))then
        write(error_unit,*) "Cannot allocate mode to '",abbrev_in,"' since it is already allocated"
        ERROR STOP
    endif
    allocate(F_mode_rank1_point::mode)
    select type(mode)
    class is(F_mode_rank1_point)
        Call mode%init_order(abbrev_in)
    end select
end subroutine
end module
