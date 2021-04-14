module m_mode_construction_rankN_eigen
use,intrinsic :: iso_c_binding
use m_mode_construction
use m_derived_types, only : lattice,number_different_order_parameters
use m_eigen_mode_interface
implicit none
private
public F_mode_rankN_eigen


type, extends(F_mode) :: F_mode_rankN_eigen
    type(C_PTR)         :: modes=c_null_ptr    !array of sparse matrices
    integer,allocatable :: ratio(:)             !ratio of dim_mat(1)/dim_mat(2) for each mode -> how many contributions per state
    contains
    !necessary routines as defined by class
    procedure   :: get_mode   !subroutine which returns the mode 
    procedure   :: get_mode_disc
    procedure   :: get_mode_exc
    procedure   :: mode_reduce_comp
    procedure   :: get_ind_site

    procedure   :: get_mode_single_cont
    procedure   :: get_mode_single_disc

    procedure   :: copy
    procedure   :: destroy
    procedure   :: is_same

    !MPI
    procedure   :: bcast
    procedure   :: send
    procedure   :: recv

    !local construction routine
    procedure   :: init_order
end type

contains


subroutine get_ind_site(this,comp,site,ind)
    class(F_mode_rankN_eigen),intent(in)   :: this
    integer,intent(in)                          :: comp  !mode index
    integer,intent(in)                          :: site    !entry
    integer,intent(inout),allocatable           :: ind(:)

    integer         :: inner_dim_mode, i

    ERROR STOP "IMPLEMENT"
end subroutine

subroutine get_mode_disc(this,lat,ind,vec)
    use, intrinsic :: iso_c_binding, only: C_PTR, C_loc
    use m_type_lattice, only: dim_modes_inner
    class(F_mode_rankN_eigen),intent(in)    :: this
    type(lattice),intent(in)                :: lat
    integer,intent(in)                      :: ind(:)
    real(8),intent(inout),allocatable       :: vec(:)

    real(8),pointer :: mode_base(:)
    type(C_PTR)     :: cptr(this%N_mode)    !pointers to modes
    integer         :: i

    if(.not.allocated(vec)) allocate(vec(size(ind)))
    do i=1,this%N_mode
        Call lat%set_order_point(this%order(i),mode_base)
        cptr(i)=C_loc(mode_base)
    enddo
    vec=1.0d0
    Call eigen_get_disc(this%modes,cptr,size(ind),ind,vec)
end subroutine

subroutine get_mode_single_disc(this,lat,comp,site,ind,vec)
    use, intrinsic :: iso_c_binding, only: C_PTR, C_loc
    use m_type_lattice, only: dim_modes_inner
    class(F_mode_rankN_eigen),intent(in)    :: this
    type(lattice),intent(in)                :: lat
    integer,intent(in)                      :: comp  !mode index
    integer,intent(in)                      :: site    !entry
    integer,intent(inout),allocatable       :: ind(:)
    real(8),intent(inout),allocatable       :: vec(:)

    integer         :: ind_in(dim_modes_inner(this%order(comp)))
    integer         :: inner_dim_mode, N_out
    real(8),pointer :: mode_base(:)
    type(C_PTR)     :: cptr(this%N_mode)    !pointers to modes
    integer         :: i

    inner_dim_mode=dim_modes_inner(this%order(comp))
    ind_in=[((site-1)*inner_dim_mode+i,i=1,inner_dim_mode)]
    N_out=inner_dim_mode*this%ratio(comp)
    if(.not.allocated(ind)) allocate(ind(N_out))
    if(.not.allocated(vec)) allocate(vec(N_out))

    do i=1,this%N_mode
        Call lat%set_order_point(this%order(i),mode_base)
        cptr(i)=C_loc(mode_base)
    enddo
    vec=1.0d0
    Call eigen_get_mode_single_disc(this%modes,cptr,comp-1,size(ind_in),ind_in,N_out,ind,vec)
    ind=ind+1
end subroutine

subroutine get_mode_single_cont(this,lat,order,i,modes,vec,bnd)
    class(F_mode_rankN_eigen),intent(in)  :: this
    type(lattice),intent(in)                    :: lat
    integer,intent(in)                          :: order
    integer,intent(in)                          :: i
    real(8),pointer,intent(out)                 :: modes(:)
    integer,intent(out)                         :: bnd(2)
    real(8),allocatable,target,intent(out)      :: vec(:)   !space to allocate array if not single operator

    ERROR STOP "IMPLEMENT"  !is this really needed?
end subroutine


subroutine get_mode_exc(this,lat,comp,vec)
    use, intrinsic :: iso_fortran_env, only : error_unit
    class(F_mode_rankN_eigen),intent(in)   :: this
    type(lattice),intent(in)                    :: lat       !lattice type which knows about all states
    integer,intent(in)                          :: comp
    real(8),intent(inout)                       :: vec(:)

    real(8)         :: tmp_internal(this%mode_size,size(this%order))
    real(8),pointer :: mode_base(:)
    integer         :: i

    tmp_internal=1.d0
    do i=1,comp-1
        Call lat%set_order_point(this%order(i),mode_base)
        Call eigen_get_mode_i(this%modes,i-1,mode_base,tmp_internal(:,i))
    enddo
    do i=comp+1,this%N_mode
        Call lat%set_order_point(this%order(i),mode_base)
        Call eigen_get_mode_i(this%modes,i-1,mode_base,tmp_internal(:,i))
    enddo
    vec=product(tmp_internal,dim=2)
    nullify(mode_base)
end subroutine

subroutine mode_reduce_comp(this,lat,vec_in,comp,vec_out)
    use, intrinsic :: iso_fortran_env, only : error_unit
    class(F_mode_rankN_eigen),intent(in)  :: this
    real(8),intent(in)                          :: vec_in(:)
    type(lattice),intent(in)                    :: lat       !lattice type which knows about all states
    integer,intent(in)                          :: comp !of which operator the first entry is kept
    real(8),intent(out)                         :: vec_out(lat%dim_modes(this%order(comp))*lat%Ncell)

    vec_out=0.0d0
    Call eigen_mode_reduce(this%modes,comp-1,vec_in,vec_out)
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
        F_out%ratio=this%ratio
        Call eigen_modes_copy(this%modes,F_out%modes) 
    class default
        ERROR STOP "FAILED TO COPY F_mode_rankN_eigen mode to F_out"
    end select
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
    Call this%init_base(order,mat(1)%dim_mat(1))
    if(size(mat)/=len(abbrev_in)) ERROR STOP "Matrix size has the be the same as the length of abbrev_in"
    allocate(this%ratio(size(mat)),source=0)

    Call eigen_alloc_arr(size(mat),this%modes)
    do i=1,size(mat)
        mat(i)%row=mat(i)%row-1
        mat(i)%col=mat(i)%col-1
        Call eigen_init_mode(this%modes,i-1,mat(i)%nnz,mat(i)%dim_mat,mat(i)%row,mat(i)%col,mat(i)%val)
        this%ratio(i)=mat(i)%dim_mat(1)/mat(i)%dim_mat(2)
        Call mat(i)%destroy()
    enddo
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!            MPI ROUTINES           !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine bcast(this,comm)
    use mpi_basic                
    class(F_mode_rankN_eigen),intent(inout) ::  this
    type(mpi_type),intent(in)               ::  comm
#ifdef CPP_MPI

    Call this%bcast_base(comm)
    ERROR STOP "IMPLEMENT"
#else
    continue
#endif
end subroutine 

subroutine send(this,ithread,tag,com)
!    use mpi_basic                
    class(F_mode_rankN_eigen),intent(in)    :: this
    integer,intent(in)          :: ithread
    integer,intent(in)          :: tag
    integer,intent(in)          :: com

#ifdef CPP_MPI
    integer     :: ierr

    Call this%send_base(ithread,tag,com)
    ERROR STOP "IMPLEMENT"
#else
    continue
#endif
end subroutine

subroutine recv(this,ithread,tag,com)
!    use mpi_basic                
    class(F_mode_rankN_eigen),intent(inout) :: this
    integer,intent(in)          :: ithread
    integer,intent(in)          :: tag
    integer,intent(in)          :: com

#ifdef CPP_MPI
    integer     :: ierr

    Call this%recv_base(ithread,tag,com)
    ERROR STOP "IMPLEMENT"
#else
    continue
#endif
end subroutine





end module
