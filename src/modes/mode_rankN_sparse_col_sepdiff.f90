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
    procedure  :: reduce_other_exc => reduce_other_exc_nosym
    procedure  :: get_mode_exc_ind
    procedure  :: mode_reduce_ind_sum
end type

contains

subroutine reduce_other_exc_nosym(this,lat,op_keep,vec_other,res)
    class(F_mode_rankN_sparse_col_sepdiff),intent(in)   :: this
    type(lattice),intent(in)                            :: lat       !lattice type which knows about all states
    integer,intent(in)                                  :: op_keep
    real(8),intent(in)                                  :: vec_other(:) !other vector to which the Hamiltonian has been multiplied
    real(8),intent(inout)                               :: res(:)

    real(8)         :: tmp(size(vec_other))   !multipied, but not reduced
    logical         :: reduce(size(this%order))
    integer         :: i

    reduce=this%order==op_keep
    res=0.0d0
    do i=1,size(this%dat)
        if(.not.reduce(i)) cycle
        Call this%get_mode_exc_ind(lat,i,tmp)
        tmp=vec_other*tmp
        Call this%mode_reduce_ind_sum(lat,tmp,i,op_keep,res)
    enddo
end subroutine

subroutine get_mode_exc_ind(this,lat,ind,vec)
    use, intrinsic :: iso_fortran_env, only : error_unit
    class(F_mode_rankN_sparse_col_sepdiff),intent(in)   :: this
    type(lattice),intent(in)                            :: lat       !lattice type which knows about all states
    integer,intent(in)                                  :: ind
    real(8),intent(inout)                               :: vec(:)

    real(8)         :: tmp_internal(this%N_mode,size(this%dat))
    real(8),pointer :: mode_base(:)
    integer         :: i

    if(size(vec)/=this%N_mode) STOP "mode exc call has wrong size for vector"
    if(ind<1.or.ind>size(this%order))then
        write(error_unit,'(//A,I6)') "Tried to get mode excluding index no.:", ind
        write(error_unit,*) "But the mode only contains the order:", this%order
        ERROR STOP "This makes no sense and should probably prevented earlier in the code"
    endif

    tmp_internal=1.d0
    do i=1,ind-1
        Call lat%set_order_point(this%order(i),mode_base)
        tmp_internal(:,i)=mode_base(this%dat(i)%col)
    enddo
    do i=ind+1,size(this%dat)
        Call lat%set_order_point(this%order(i),mode_base)
        tmp_internal(:,i)=mode_base(this%dat(i)%col)
    enddo
    vec=product(tmp_internal,dim=2)
    nullify(mode_base)
end subroutine

subroutine mode_reduce_ind_sum(this,lat,vec_in,ind,opind,vec_out)
    use, intrinsic :: iso_fortran_env, only : error_unit
    class(F_mode_rankN_sparse_col_sepdiff),intent(in)   :: this
    real(8),intent(in)                          :: vec_in(:)
    type(lattice),intent(in)                    :: lat       !lattice type which knows about all states
    integer,intent(in)                          :: ind
    integer,intent(in)                          :: opind
    real(8),intent(inout)                       :: vec_out(lat%dim_modes(opind)*lat%Ncell)

    integer     ::  i

    do i=1,this%dat(ind)%nnz
        vec_out(this%dat(ind)%col(i))=vec_out(this%dat(ind)%col(i))+vec_in(i)
    enddo
end subroutine
end module
