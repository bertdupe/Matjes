module m_mode_construction
! module that contains functions which combines the different sites on the left and the right of the Hamiltonians
! to transform Hamiltonians of ranks >2 to Hamiltonians of rank 2
use m_derived_types, only : lattice
implicit none

type, abstract :: F_mode
    contains
    procedure(int_get_mode),deferred          :: get_mode   !subroutine which returns the mode 
    !TODO: something which compares the mode construction when checking of Hamiltonian terms can be combined
end type

abstract interface
    subroutine int_get_mode(this,lat,mode,tmp)
        import F_mode,lattice
        class(F_mode),intent(in)                   :: this
        type(lattice),intent(in)                   :: lat       !lattice type which knows about all states
        real(8),intent(out),pointer                :: mode(:)   !pointer to required mode
        real(8),allocatable,target,intent(inout)   :: tmp(:)    !possible temporary storage for mode pointer
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
    procedure   :: init_order
    procedure   :: get_mode   !subroutine which returns the mode 
end type

contains

subroutine get_mode(this,lat,mode,tmp)
    class(F_mode_rank1_point),intent(in)       :: this
    type(lattice),intent(in)                   :: lat       !lattice type which knows about all states
    real(8),intent(out),pointer                :: mode(:)   !pointer to required mode
    real(8),allocatable,target,intent(inout)   :: tmp(:)    !not used here

    Call lat%set_order_point(this%order,mode)
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

