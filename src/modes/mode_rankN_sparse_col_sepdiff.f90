module m_mode_construction_rankN_sparse_col_sepdiff
!make the derivative for each entry of the derived order separately
use m_mode_construction
use m_derived_types, only : lattice,number_different_order_parameters
use m_coo_mat
use m_mode_construction_rankN_sparse_col
implicit none
private
public F_mode_rankN_sparse_col_sepdiff

type, extends(F_mode_rankN_sparse_col) :: F_mode_rankN_sparse_col_sepdiff
    contains
    procedure  :: reduce_other_exc
    procedure  :: mode_reduce_ind_sum
end type

contains

subroutine reduce_other_exc(this,lat,op_keep,vec_other,res)
    class(F_mode_rankN_sparse_col_sepdiff),intent(in)   :: this
    type(lattice),intent(in)                            :: lat       !lattice type which knows about all states
    integer,intent(in)                                  :: op_keep
    real(8),intent(in)                                  :: vec_other(:) !other vector to which the Hamiltonian has been multiplied
    real(8),intent(inout)                               :: res(:)

    real(8)         :: tmp(size(vec_other))   !multipied, but not reduced
    logical         :: reduce(this%N_mode)
    integer         :: i

    reduce=this%order==op_keep
    res=0.0d0
    do i=1,size(this%dat)
        if(.not.reduce(i)) cycle
        Call this%get_mode_exc_ind(lat,i,tmp)
        tmp=vec_other*tmp
        Call this%mode_reduce_ind_sum(lat,tmp,i,res)
    enddo
end subroutine

subroutine mode_reduce_ind_sum(this,lat,vec_in,ind,vec_out)
    !like mode_reduce_ind, but sums vec_out up
    use, intrinsic :: iso_fortran_env, only : error_unit
    class(F_mode_rankN_sparse_col_sepdiff),intent(in)   :: this
    real(8),intent(in)          :: vec_in(:)
    type(lattice),intent(in)    :: lat       !lattice type which knows about all states
    integer,intent(in)          :: ind
    real(8),intent(inout)       :: vec_out(lat%dim_modes(this%order(ind))*lat%Ncell)
    integer     ::  i
    do i=1,this%dat(ind)%nnz
        vec_out(this%dat(ind)%col(i))=vec_out(this%dat(ind)%col(i))+vec_in(i)
    enddo
end subroutine
end module
