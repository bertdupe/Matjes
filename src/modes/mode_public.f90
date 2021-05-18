module m_mode_public
use m_derived_types,only: lattice
use m_mode_construction, only: F_mode
use m_mode_construction_rank1_point, only: F_mode_rank1_point
use m_mode_construction_rankN_full_manual, only: F_mode_rankN_full_manual
#ifdef CPP_EIGEN
use m_mode_construction_rankN_eigen, only: F_mode_rankN_eigen
#endif
use m_mode_construction_rankN_sparse_col, only: F_mode_rankN_sparse_col
use m_mode_construction_rankN_sparse_col_sepdiff, only: F_mode_rankN_sparse_col_sepdiff
#if 0
use m_mode_construction_rankN_sparse_coo, only: F_mode_rankN_sparse_coo
use m_mode_construction_rankN_sparse_coo_red, only: F_mode_rankN_sparse_coo_red
#endif
public

contains 

subroutine get_mode_ident(mode,ident)
    use, intrinsic  ::  ISO_FORTRAN_ENV, only: error_unit
    !subroutine which return an identifier to allow the allocation of the correct type encoded in an integer
    class(F_mode),allocatable,intent(in)    :: mode
    integer,intent(out)     ::  ident

    if(.not.allocated(mode))then
        write(error_unit,"(///A)") "Cannot get mode indentifier, as it is not allocated"
        ERROR STOP
    endif
    select type(mode)
    type is(F_mode_rank1_point)
        ident=1
    type is(F_mode_rankN_full_manual)
        ident=2
    type is(F_mode_rankN_sparse_col)
        ident=3
    type is(F_mode_rankN_sparse_col_sepdiff)
        ident=4
#ifdef CPP_EIGEN
    type is(F_mode_rankN_eigen)
        ident=5
#else
    ERROR STOP "CANNOT USE EIGEN MODE DESCIPTION WITHOUT EIGEN LIBRARY"
#endif
    class default
        ident=0
        write(error_unit,"(///A)") "Failed to identify type of mode, probably a programming mistake"
        ERROR STOP
    end select
end subroutine

subroutine set_mode_ident(mode,ident)
    use, intrinsic  ::  ISO_FORTRAN_ENV, only: error_unit
    !subroutine which allocates the mode to the correct type encoded by the integer of from get_mode_ident
    class(F_mode),allocatable,intent(inout) :: mode
    integer,intent(in)      ::  ident

    if(allocated(mode))then
        write(error_unit,"(///A)") "Called set_mode_ident with already allocated mode which should not"
        ERROR STOP
    endif
    select case(ident)
    case(1)
        allocate(F_mode_rank1_point::mode)
    case(2)
        allocate(F_mode_rankN_full_manual::mode)
    case(3)
        allocate(F_mode_rankN_sparse_col::mode)
    case(4)
        allocate(F_mode_rankN_sparse_col_sepdiff::mode)
#ifdef CPP_EIGEN
    case(5)
        allocate(F_mode_rankN_eigen::mode)
#else
    ERROR STOP "CANNOT USE EIGEN MODE DESCIPTION WITHOUT EIGEN LIBRARY"
#endif
    case default
        write(error_unit,"(///A,I6,/,A)") "Received unexpected identified to allocate mode: ",ident,"Was an mode implementation forgotten?"
        ERROR STOP
    end select
end subroutine



subroutine mode_set_rank1(mode,lat,abbrev_in)
    !function to set mode to single rank avoiding ecasing all the allocation and type select stuff to be more concise in the Energy definitions
    use, intrinsic :: iso_fortran_env, only : error_unit
    class(F_mode),allocatable,intent(inout) :: mode
    type(lattice),intent(in)                :: lat       !lattice type which knows about all states
    character(len=1), intent(in)            :: abbrev_in

    if(allocated(mode))then
        write(error_unit,*) "Cannot allocate mode to '",abbrev_in,"' since it is already allocated"
        ERROR STOP
    endif
    allocate(F_mode_rank1_point::mode)
    select type(mode)
    class is(F_mode_rank1_point)
        Call mode%init_order(lat,abbrev_in)
    end select
end subroutine

subroutine mode_set_rankN(mode,abbrev_in,lat,implementation)
    !might not work after all, since many different input parameters might be necessary
    !function to set mode to single rank avoiding ecasing all the allocation and type select stuff to be more concise in the Energy definitions
    use, intrinsic :: iso_fortran_env, only : error_unit
    class(F_mode),allocatable,intent(inout) :: mode
    type(lattice),intent(in)                :: lat       !lattice type which knows about all states
    character(len=*), intent(in)            :: abbrev_in
    integer                                 :: implementation   !integer input to choose between different implementations 

    if(allocated(mode))then
        write(error_unit,*) "Cannot allocate mode to '",abbrev_in,"' since it is already allocated"
        ERROR STOP
    endif
    if(implementation==1)then
        allocate(F_mode_rankN_full_manual::mode)
        select type(mode)
        type is(F_mode_rankN_full_manual)
            Call mode%init_order(lat,abbrev_in)
        end select
    else
        write(error_unit,*) "Trying to set rankN mode=",implementation
        write(error_unit,*) "This is not implemented and probably a programming mistake"
        ERROR STOP
    endif
end subroutine

subroutine mode_set_rankN_sparse(mode,abbrev_in,lat,mat,implementation)
    !might not work after all, since many different input parameters might be necessary
    !function to set mode to single rank avoiding ecasing all the allocation and type select stuff to be more concise in the Energy definitions
    use, intrinsic :: iso_fortran_env, only : error_unit
    use m_coo_mat, only: coo_mat
    class(F_mode),allocatable,intent(inout) :: mode
    type(lattice),intent(in)                :: lat       !lattice type which knows about all states
    character(len=*), intent(in)            :: abbrev_in
    type(coo_mat),intent(inout)             :: mat(:)    !input matrices, destroyed when returned

    integer                                 :: implementation   !integer input to choose between different implementations 

    if(allocated(mode))then
        write(error_unit,*) "Cannot allocate mode to '",abbrev_in,"' since it is already allocated"
        ERROR STOP
    endif
    if(implementation==1)then
        !save the mode construction procedure only by column array( assume val=1 & rows consecutive)
        allocate(F_mode_rankN_sparse_col::mode)
        select type(mode)
        type is(F_mode_rankN_sparse_col)
            Call mode%init_order(lat,abbrev_in,mat)
        end select
    elseif(implementation==2)then
        !save the mode construction procedure only by column array( assume val=1 & rows consecutive)
        !does the derivative operation for each mode with the operator separately 
        allocate(F_mode_rankN_sparse_col_sepdiff::mode)
        select type(mode)
        type is(F_mode_rankN_sparse_col_sepdiff)
            Call mode%init_order(lat,abbrev_in,mat)
        end select
    elseif(implementation==3)then
#ifdef CPP_EIGEN
        !slower than using *_col type
        !saves the construction procedure in a eigen sparse matrix 
        allocate(F_mode_rankN_eigen::mode)
        select type(mode)
        type is(F_mode_rankN_eigen)
            Call mode%init_order(lat,abbrev_in,mat)
        end select
#else
    ERROR STOP "CANNOT USE EIGEN MODE DESCIPTION WITHOUT EIGEN LIBRARY"
#endif
#if 0
    elseif(implementation==98)then
        !rather obsolete, unless *_col does not work
        !save the mode construction procedure as coo_mat type with col,row, val
        allocate(F_mode_rankN_sparse_coo::mode)
        select type(mode)
        type is(F_mode_rankN_sparse_coo)
            Call mode%init_order(lat,abbrev_in,mat)
        end select
    elseif(implementation==99)then
        !rather obsolete, unless *_col does not work
        !save the mode construction procedure as coo_mat type with col,row, val
        !uses only col for construction (assume val=1 & rows consecutive)
        allocate(F_mode_rankN_sparse_coo_red::mode)
        select type(mode)
        type is(F_mode_rankN_sparse_coo_red)
            Call mode%init_order(lat,abbrev_in,mat)
        end select
#endif
    else
        write(error_unit,*) "Trying to set rankN mode=",implementation
        write(error_unit,*) "This is not implemented and probably a programming mistake"
        ERROR STOP
    endif
end subroutine

end module
