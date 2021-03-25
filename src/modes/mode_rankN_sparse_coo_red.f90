module m_mode_construction_rankN_sparse_coo_red
!assume all values are 1 and row is just [(i,i=1,Nmode)]
use m_mode_construction
use m_derived_types, only : lattice,number_different_order_parameters
use m_coo_mat
use m_mode_construction_rankN_sparse_coo
implicit none
private
public F_mode_rankN_sparse_coo_red

type, extends(F_mode_rankN_sparse_coo) :: F_mode_rankN_sparse_coo_red !contains all entries
    contains
    !necessary routines as defined by class
    procedure   :: get_mode   !subroutine which returns the mode 
    procedure   :: get_mode_exc_ind
    procedure   :: mode_reduce_ind 
end type

contains

subroutine get_mode_exc_ind(this,lat,ind,vec)
    use, intrinsic :: iso_fortran_env, only : error_unit
    class(F_mode_rankN_sparse_coo_red),intent(in)   :: this
    type(lattice),intent(in)                        :: lat       !lattice type which knows about all states
    integer,intent(in)                              :: ind
    real(8),intent(inout)                           :: vec(:)

    real(8)         :: tmp_internal(this%mode_size,size(this%mat))
    logical         :: exclude(size(this%order))
    real(8),pointer :: mode_base(:)
    integer         :: i,j

    tmp_internal=1.d0
    do i=1,ind-1
        Call lat%set_order_point(this%order(i),mode_base)
        tmp_internal(:,i)=mode_base(this%mat(i)%col)
    enddo
    do i=ind+1,this%N_mode
        Call lat%set_order_point(this%order(i),mode_base)
        tmp_internal(:,i)=mode_base(this%mat(i)%col)
    enddo
    vec=product(tmp_internal,dim=2)
    nullify(mode_base)
end subroutine

subroutine mode_reduce_ind(this,lat,vec_in,ind,vec_out)
    use, intrinsic :: iso_fortran_env, only : error_unit
    class(F_mode_rankN_sparse_coo_red),intent(in)   :: this
    real(8),intent(in)                              :: vec_in(:)
    type(lattice),intent(in)                        :: lat       !lattice type which knows about all states
    integer,intent(in)                              :: ind
    real(8),intent(out)                             :: vec_out(lat%dim_modes(this%order(ind))*lat%Ncell)

    integer     ::  i
    vec_out=0.0d0
    do i=1,this%mat(ind)%nnz
        vec_out(this%mat(ind)%col(i))=vec_out(this%mat(ind)%col(i))+vec_in(i)
    enddo
end subroutine

subroutine get_mode(this,lat,mode,tmp)
    class(F_mode_rankN_sparse_coo_red),intent(in)   :: this
    type(lattice),intent(in)                        :: lat       !lattice type which knows about all states
    real(8),intent(out),pointer                     :: mode(:)   !pointer to required mode
    real(8),allocatable,target,intent(inout)        :: tmp(:)

    integer         :: N
    real(8)         :: tmp_internal(this%mode_size,size(this%mat))
    real(8),pointer :: mode_base(:)

    integer         :: i,j

    allocate(tmp(this%mode_size),source=0.0d0)
    mode=>tmp

    tmp_internal=0.d0
    do i=1,size(this%mat)
        Call lat%set_order_point(this%order(i),mode_base)
        tmp_internal(:,i)=mode_base(this%mat(i)%col)
    enddo
    mode=product(tmp_internal,dim=2)
    nullify(mode_base)
end subroutine
end module
