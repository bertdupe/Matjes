module m_mode_construction_rankN_eigen
use,intrinsic :: iso_c_binding
use m_mode_construction
use m_derived_types, only : lattice,number_different_order_parameters
use m_eigen_mode_interface
implicit none
private
public F_mode_rankN_eigen


type, extends(F_mode) :: F_mode_rankN_eigen
    type(C_PTR)     :: modes=c_null_ptr    !array of sparse matrices
    integer         :: mode_size=-1
    contains
    !necessary routines as defined by class
    procedure   :: get_mode   !subroutine which returns the mode 
    procedure   :: get_mode_exc_ind
    procedure   :: mode_reduce_ind
    procedure   :: get_mode_single_cont

    procedure   :: copy
    procedure   :: bcast
    procedure   :: destroy
    procedure   :: is_same
    !local construction routine
    procedure   :: init_order
end type

contains


subroutine get_mode_single_disc(this,lat,ind_mode,site,ind,vec)
    class(F_mode_rankN_eigen),intent(in)    :: this
    type(lattice),intent(in)                :: lat
    integer,intent(in)                      :: ind_mode  !mode index
    integer,intent(in)                      :: site    !entry
    integer,intent(inout),allocatable       :: ind(:)
    real(8),intent(inout),allocatable       :: vec(:)

    real(8),pointer :: mode_base(:)
!    real(8)         :: tmp_internal(N,size(this%order))

!    Call eigen_mode_get_ind
    ERROR STOP "IMPLEMENT"

end subroutine

subroutine get_mode_single_cont(this,lat,order,i,modes,vec,bnd)
    class(F_mode_rankN_eigen),intent(in)  :: this
    type(lattice),intent(in)                    :: lat
    integer,intent(in)                          :: order
    integer,intent(in)                          :: i
    real(8),pointer,intent(out)                 :: modes(:)
    integer,intent(out)                         :: bnd(2)
    real(8),allocatable,target,intent(out)      :: vec(:)   !space to allocate array if not single operator

    ERROR STOP "IMPLEMENT"
end subroutine


subroutine get_mode_exc_ind(this,lat,ind,vec)
    use, intrinsic :: iso_fortran_env, only : error_unit
    class(F_mode_rankN_eigen),intent(in)   :: this
    type(lattice),intent(in)                    :: lat       !lattice type which knows about all states
    integer,intent(in)                          :: ind
    real(8),intent(inout)                       :: vec(:)

    real(8)         :: tmp_internal(this%mode_size,size(this%order))
    real(8),pointer :: mode_base(:)
    integer         :: i

    tmp_internal=1.d0
    do i=1,ind-1
        Call lat%set_order_point(this%order(i),mode_base)
        Call eigen_get_mode_i(this%modes,i-1,mode_base,tmp_internal(:,i))
    enddo
    do i=ind+1,this%N_mode
        Call lat%set_order_point(this%order(i),mode_base)
        Call eigen_get_mode_i(this%modes,i-1,mode_base,tmp_internal(:,i))
    enddo
    vec=product(tmp_internal,dim=2)
    nullify(mode_base)
end subroutine

subroutine mode_reduce_ind(this,lat,vec_in,ind,vec_out)
    use, intrinsic :: iso_fortran_env, only : error_unit
    class(F_mode_rankN_eigen),intent(in)  :: this
    real(8),intent(in)                          :: vec_in(:)
    type(lattice),intent(in)                    :: lat       !lattice type which knows about all states
    integer,intent(in)                          :: ind !of which operator the first entry is kept
    real(8),intent(out)                         :: vec_out(lat%dim_modes(this%order(ind))*lat%Ncell)

    vec_out=0.0d0
    Call eigen_mode_reduce(this%modes,ind-1,vec_in,vec_out)
end subroutine


subroutine get_mode(this,lat,mode,tmp)
    class(F_mode_rankN_eigen),intent(in)   :: this
    type(lattice),intent(in)                    :: lat       !lattice type which knows about all states
    real(8),intent(out),pointer                 :: mode(:)   !pointer to required mode
    real(8),allocatable,target,intent(inout)    :: tmp(:)

    real(8)         :: tmp_internal(this%mode_size,size(this%order))
    real(8),pointer :: mode_base(:)
    integer         :: i

    allocate(tmp(this%mode_size),source=0.0d0)
    mode=>tmp
!
    do i=1,size(this%order)
        Call lat%set_order_point(this%order(i),mode_base)
        Call eigen_get_mode_i(this%modes,i-1,mode_base,tmp_internal(:,i))
    enddo
    mode=product(tmp_internal,dim=2)
    nullify(mode_base)
end subroutine

function is_same(this,comp)result(same)
    class(F_mode_rankN_eigen),intent(in)  :: this
    class(F_mode),intent(in)              :: comp
    logical                               :: same

    ERROR STOP "IMPLEMENT"
end function

subroutine destroy(this)
    class(F_mode_rankN_eigen),intent(inout) ::  this
    integer ::  i

    Call eigen_mode_destroy(this%modes)
end subroutine

subroutine copy(this,F_out)
    class(F_mode_rankN_eigen),intent(in)    :: this
    class(F_mode),allocatable,intent(inout) :: F_out

    Call this%copy_base(F_out)
    select type(F_out)
    class is(F_mode_rankN_eigen)
        F_out%mode_size=this%mode_size
        Call eigen_modes_copy(this%modes,F_out%modes) 
    class default
        ERROR STOP "FAILED TO COPY F_mode_rankN_eigen mode to F_out"
    end select
end subroutine

subroutine bcast(this,comm)
    use mpi_basic                
    class(F_mode_rankN_eigen),intent(inout) ::  this        !this might fail if the server threads non-allocated class(F_mode), TAKE CARE OF THIS IN HAM_BASE
    type(mpi_type),intent(in)               ::  comm
#ifdef CPP_MPI

    Call bcast_base(this,comm)
    ERROR STOP "IMPLEMENT"
#else
    continue
#endif
end subroutine 

subroutine init_order(this,lat,abbrev_in,mat)
    use m_derived_types, only: op_abbrev_to_int
    use m_coo_mat, only: coo_mat
    class(F_mode_rankN_eigen),intent(inout) :: this
    type(lattice),intent(in)                :: lat       !lattice type which knows about all states
    character(len=*), intent(in)            :: abbrev_in !considered order abbreviations
    type(coo_mat),intent(inout)             :: mat(:)    !input matrices, destroyed when returned
    integer     :: order(len(abbrev_in))
    integer     :: i     
    integer     :: order_occ(number_different_order_parameters)


    order=op_abbrev_to_int(abbrev_in)
    Call this%init_base(order)
    if(size(mat)/=len(abbrev_in)) ERROR STOP "Matrix size has the be the same as the length of abbrev_in"
    this%mode_size=mat(1)%dim_mat(1)

    Call eigen_alloc_arr(size(mat),this%modes)
    do i=1,size(mat)
        mat(i)%row=mat(i)%row-1
        mat(i)%col=mat(i)%col-1
        Call eigen_init_mode(this%modes,i-1,mat(i)%nnz,mat(i)%dim_mat,mat(i)%row,mat(i)%col,mat(i)%val)
        Call mat(i)%destroy()
    enddo
end subroutine
end module
