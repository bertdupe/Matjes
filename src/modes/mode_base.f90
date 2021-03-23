module m_mode_construction
! module that contains functions which combines the different sites on the left and the right of the Hamiltonians
! to transform Hamiltonians of ranks >2 to Hamiltonians of rank 2
use m_derived_types, only : lattice, number_different_order_parameters
implicit none

private
public F_mode

type, abstract :: F_mode
    integer             :: order_occ(number_different_order_parameters)=0
    contains
    procedure(int_get_mode),deferred        :: get_mode       !subroutine which returns the mode 
    procedure(int_get_mode_exc),deferred    :: get_mode_exc   !subroutine which returns the mode excluding one order parameter 
    procedure(int_mode_reduce),deferred     :: mode_reduce    !subroutine which reduces an input vector to the shape of a single constituent state according the the F_mode rules

    !procedure(int_get_mode_single),deferred :: get_mode_single       !subroutine which returns the mode 
    !procedure(int_get_mode_exc),deferred    :: get_mode_exc   !subroutine which returns the mode excluding one order parameter 
    !procedure(int_mode_reduce),deferred     :: mode_reduce    !subroutine which reduces an input vector to the shape of a single constituent state according the the F_mode rules

    procedure(int_is_same),deferred         :: is_same
    procedure(int_destroy),deferred         :: destroy
    procedure(int_copy),deferred            :: copy
    procedure(int_bcast),deferred           :: bcast

    procedure                               :: init_base
    procedure                               :: copy_base

    procedure                               :: reduce_other_exc => reduce_other_exc_default
end type

abstract interface
    subroutine int_get_mode(this,lat,mode,tmp)
        !subroutine which return a pointer to the mode, either saved in the tmp array or at some other independent array (depending on the implementation)
        import F_mode,lattice
        class(F_mode),intent(in)                   :: this
        type(lattice),intent(in)                   :: lat       !lattice type which knows about all states
        real(8),intent(out),pointer                :: mode(:)   !pointer to required mode
        real(8),allocatable,target,intent(inout)   :: tmp(:)    !possible temporary storage for mode pointer
    end subroutine

    subroutine int_get_mode_exc(this,lat,op_exc,vec)
        !subroutine which gets the mode, excluding the first occation of the state indexed by op_exc
        import F_mode,lattice
        class(F_mode),intent(in)                    :: this
        type(lattice),intent(in)                    :: lat       !lattice type which knows about all states
        integer,intent(in)                          :: op_exc    !operator index which is not multiplied [1,number_different_order_parameters]
        real(8),intent(inout)                       :: vec(:) !result mode excluding the first state with op_exc
    end subroutine

    subroutine int_mode_reduce(this,lat,vec_in,op_keep,vec_out)
        !subroutine which reduces the input mode according to the rules of this F_modes construction rules
        ! to be in the basis of the first mode with the op_keep index
        import F_mode,lattice
        class(F_mode),intent(in)        :: this
        real(8),intent(in)              :: vec_in(:)
        type(lattice),intent(in)        :: lat       !lattice type which knows about all states
        integer,intent(in)              :: op_keep   !of which operator the first entry is kept
        real(8),intent(out)             :: vec_out(lat%dim_modes(op_keep)*lat%Ncell)
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


subroutine init_base(this,order_occ)
        class(F_mode),intent(inout) :: this
        integer,intent(in)          :: order_occ(number_different_order_parameters)
        this%order_occ=order_occ
end subroutine

subroutine copy_base(this,F_out)
        class(F_mode),intent(in)                :: this
        class(F_mode),intent(inout),allocatable :: F_out
        if(.not.allocated(F_out)) allocate(F_out,mold=this)
        F_out%order_occ=this%order_occ
end subroutine

subroutine reduce_other_exc_default(this,lat,op_keep,vec_other,res)
    !subroutine to get the derivative, may be overwritten
    class(F_mode),intent(in)        :: this
    type(lattice),intent(in)        :: lat       !lattice type which knows about all states
    integer,intent(in)              :: op_keep  !operator index which is to be kept (1,number_different_order_parameters)
    real(8),intent(in)              :: vec_other(:) !other vector to which the Hamiltonian has been multiplied
    real(8),intent(inout)           :: res(:)

    real(8)         :: tmp(size(vec_other))
    integer         :: i

    Call this%get_mode_exc(lat,op_keep,tmp)
    tmp=vec_other*tmp
    Call this%mode_reduce(lat,tmp,op_keep,res)
    tmp=tmp*real(this%order_occ(op_keep),8)
end subroutine

end module 
