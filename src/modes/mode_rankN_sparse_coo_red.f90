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
    procedure   :: get_mode_exc
    procedure   :: mode_reduce  
end type

contains

subroutine get_mode_exc(this,lat,op_exc,vec)
    use, intrinsic :: iso_fortran_env, only : error_unit
    class(F_mode_rankN_sparse_coo_red),intent(in)   :: this
    type(lattice),intent(in)                        :: lat       !lattice type which knows about all states
    integer,intent(in)                              :: op_exc !of which operator the first entry is kept
    real(8),intent(inout)                           :: vec(:)

    real(8)         :: tmp_internal(this%N_mode,size(this%mat))
    logical         :: exclude(size(this%order))
    real(8),pointer :: mode_base(:)
    integer         :: i,j

    if(size(vec)/=this%N_mode) STOP "mode exc call has wrong size for vector"
    exclude=.false.
    i=findloc(this%order,op_exc,dim=1)
    if(i<1.or.i>size(this%order))then
        write(error_unit,'(//A,I6)') "Tried to get mode excluding order no.:", op_exc
        write(error_unit,*) "But the mode only contains the order:", this%order
        ERROR STOP "This makes no sense and should probably prevented earlier in the code"
    endif
    exclude(i)=.true.

    tmp_internal=1.d0
    do i=1,size(this%mat)
        if(.not.exclude(i))then
            Call lat%set_order_point(this%order(i),mode_base)
            tmp_internal(:,i)=mode_base(this%mat(i)%col)
        endif
    enddo
    vec=product(tmp_internal,dim=2)
    nullify(mode_base)
end subroutine

subroutine mode_reduce(this,lat,vec_in,op_keep,vec_out)
    use, intrinsic :: iso_fortran_env, only : error_unit
    class(F_mode_rankN_sparse_coo_red),intent(in)   :: this
    real(8),intent(in)                              :: vec_in(:)
    type(lattice),intent(in)                        :: lat       !lattice type which knows about all states
    integer,intent(in)                              :: op_keep   !of which operator the first entry is kept
    real(8),intent(out)                             :: vec_out(lat%dim_modes(op_keep)*lat%Ncell)

    integer     ::  i_order,i

    i_order=findloc(this%order,op_keep,dim=1)
    vec_out=0.0d0
    do i=1,this%mat(i_order)%nnz
        vec_out(this%mat(i_order)%col(i))=vec_out(this%mat(i_order)%col(i))+vec_in(i)
    enddo
end subroutine


subroutine get_mode(this,lat,mode,tmp)
    class(F_mode_rankN_sparse_coo_red),intent(in)   :: this
    type(lattice),intent(in)                        :: lat       !lattice type which knows about all states
    real(8),intent(out),pointer                     :: mode(:)   !pointer to required mode
    real(8),allocatable,target,intent(inout)        :: tmp(:)

    integer         :: N
    real(8)         :: tmp_internal(this%N_mode,size(this%mat))
    real(8),pointer :: mode_base(:)

    integer         :: i,j

    allocate(tmp(this%N_mode),source=0.0d0)
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
