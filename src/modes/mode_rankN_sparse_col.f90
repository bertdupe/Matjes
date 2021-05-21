module m_mode_construction_rankN_sparse_col
use m_mode_construction
use m_type_lattice, only: dim_modes_inner,  lattice,number_different_order_parameters
use m_coo_mat
implicit none
private
public F_mode_rankN_sparse_col, col_mat

type col_mat
    integer :: dim_mat(2)=0
    integer :: nnz=0
    integer,allocatable :: col(:)
    integer,allocatable :: reverse(:,:)
contains
    procedure   :: bcast=> col_mat_bcast
    procedure   :: recv => col_mat_recv
    procedure   :: send => col_mat_send
end type

type, extends(F_mode) :: F_mode_rankN_sparse_col 
    type(col_mat),allocatable   :: dat(:)
    integer,allocatable         :: ratio(:)             !ratio of dim_mat(1)/dim_mat(2) for each mode -> how many contributions per state
    contains
    !necessary routines as defined by class
    procedure   :: get_mode   !subroutine which returns the mode 
    procedure   :: get_mode_disc
    procedure   :: get_mode_disc_expl
    procedure   :: get_mode_exc
    procedure   :: get_mode_exc_disc
    procedure   :: mode_reduce_comp
    procedure   :: mode_reduce_comp_disc

    procedure   :: reduce_site_vec
    procedure   :: get_mode_single_size

    procedure   :: get_ind_site
    procedure   :: get_ind_site_expl

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

subroutine reduce_site_vec(this,comp,vec_in,vec_out)
    !reduce the vector along the comp mode if the vec_in has been set with get_ind_site
    class(F_mode_rankN_sparse_col),intent(in)   :: this
    integer,intent(in)                          :: comp
    real(8),intent(in)                          :: vec_in(:)
    real(8),intent(out)                         :: vec_out(:)

#ifdef CPP_DEBUG
    if(comp>this%N_mode) ERROR STOP "COMP shall not be larger than N_mode"
#endif
    Call reduce(this%ratio(comp),dim_modes_inner(this%order(comp)),vec_in,vec_out)
end subroutine

subroutine reduce(N1,N2,vec_in,vec_out)
    integer,intent(in)      :: N1
    integer,intent(in)      :: N2
    real(8),intent(in)      :: vec_in(N1,N2)
    real(8),intent(out)     :: vec_out(N2)

    vec_out=sum(vec_in,1)
end subroutine

subroutine get_mode_single_size(this,order,dim_mode)
    !returns the size of the vector necessary to get to mode set by a single site
    class(F_mode_rankN_sparse_col),intent(in)    :: this
    integer,intent(in)          :: order
    integer,intent(out)         :: dim_mode

    integer :: order_ind,i

    if(.not.any(order/=this%order))then
#ifdef CPP_DEBUG
        write(error_unit,'(A)') "trying to get single mode of order which is not the order of the F_mode"
#endif
        dim_mode=0
    else
        order_ind=findloc(order==this%order,.true.,dim=1)
        dim_mode=dim_modes_inner(order)*this%ratio(order_ind)
    endif
end subroutine


subroutine get_ind_site(this,comp,site,ind)
    class(F_mode_rankN_sparse_col),intent(in)   :: this
    integer,intent(in)                          :: comp  !mode index
    integer,intent(in)                          :: site    !entry
    integer,intent(inout),allocatable           :: ind(:)

    integer         :: ind_site(dim_modes_inner(this%order(comp)))
    integer         :: inner_dim_mode, ratio, N_out, i

    ratio=this%ratio(comp)
    inner_dim_mode=dim_modes_inner(this%order(comp))
    ind_site=[((site-1)*inner_dim_mode+i,i=1,inner_dim_mode)]
    N_out=inner_dim_mode*ratio
    if(.not.allocated(ind)) allocate(ind(N_out),source=0)

    do i=1,inner_dim_mode
        ind((i-1)*ratio+1:i*ratio)=this%dat(comp)%reverse(:,ind_site(i))
    enddo
end subroutine

subroutine get_ind_site_expl(this,comp,site,size_out,ind)
    !get the indices corresponding to the 
    class(F_mode_rankN_sparse_col),intent(in)   :: this
    integer,intent(in)                          :: comp
    integer,intent(in)                          :: site    !entry
    integer,intent(in)                          :: size_out
    integer,intent(out)                         :: ind(size_out)

    integer         :: inner_dim_mode, ratio, i, offset

    ratio=this%ratio(comp)
    inner_dim_mode=dim_modes_inner(this%order(comp))
    offset=(site-1)*inner_dim_mode
    do i=1,inner_dim_mode
        ind((i-1)*ratio+1:i*ratio)=this%dat(comp)%reverse(:,offset+i)
    enddo
end subroutine

subroutine get_mode_disc_expl(this,lat,N,ind,vec)
    use, intrinsic :: iso_c_binding, only: C_PTR, C_loc
    use m_type_lattice, only: dim_modes_inner
    class(F_mode_rankN_sparse_col),intent(in)    :: this
    type(lattice),intent(in)                :: lat
    integer,intent(in)                      :: N
    integer,intent(in)                      :: ind(N)
    real(8),intent(out)                     :: vec(N)

    real(8),pointer :: mode_base(:)
    integer ::  i

    vec=1.0
    do i=1,size(this%dat)
        Call lat%set_order_point(this%order(i),mode_base)
        vec=vec*mode_base(this%dat(i)%col(ind))
    enddo
end subroutine

subroutine get_mode_disc(this,lat,ind,vec)
    use, intrinsic :: iso_c_binding, only: C_PTR, C_loc
    use m_type_lattice, only: dim_modes_inner
    class(F_mode_rankN_sparse_col),intent(in)    :: this
    type(lattice),intent(in)                :: lat
    integer,intent(in)                      :: ind(:)
    real(8),intent(inout),allocatable       :: vec(:)

    real(8)         :: tmp_internal(size(ind),size(this%dat))
    real(8),pointer :: mode_base(:)
    integer ::  i

    if(.not.allocated(vec)) allocate(vec(size(ind)))
    tmp_internal=1.d0
    do i=1,size(this%dat)
        Call lat%set_order_point(this%order(i),mode_base)
        tmp_internal(:,i)=mode_base(this%dat(i)%col(ind))
    enddo
    vec=product(tmp_internal,dim=2)
    nullify(mode_base)
end subroutine

subroutine get_mode_exc_disc(this,lat,comp,ind,vec)
    use, intrinsic :: iso_fortran_env, only : error_unit
    class(F_mode_rankN_sparse_col),intent(in)   :: this
    type(lattice),intent(in)                    :: lat       !lattice type which knows about all states
    integer,intent(in)                          :: comp
    integer,intent(in)                          :: ind(:)
    real(8),intent(inout)                       :: vec(:)

    real(8)         :: tmp_internal(size(vec),this%N_mode)
    real(8),pointer :: mode_base(:)
    integer         :: i

    tmp_internal=1.d0
    do i=1,comp-1
        Call lat%set_order_point(this%order(i),mode_base)
        tmp_internal(:,i)=mode_base(this%dat(i)%col(ind))
    enddo
    do i=comp+1,this%N_mode
        Call lat%set_order_point(this%order(i),mode_base)
        tmp_internal(:,i)=mode_base(this%dat(i)%col(ind))
    enddo
    vec=product(tmp_internal,dim=2)
    nullify(mode_base)
end subroutine

subroutine get_mode_exc(this,lat,comp,vec)
    use, intrinsic :: iso_fortran_env, only : error_unit
    class(F_mode_rankN_sparse_col),intent(in)   :: this
    type(lattice),intent(in)                    :: lat       !lattice type which knows about all states
    integer,intent(in)                          :: comp
    real(8),intent(inout)                       :: vec(:)

    real(8)         :: tmp_internal(this%mode_size,this%N_mode)
    real(8),pointer :: mode_base(:)
    integer         :: i

    tmp_internal=1.d0
    do i=1,comp-1
        Call lat%set_order_point(this%order(i),mode_base)
        tmp_internal(:,i)=mode_base(this%dat(i)%col)
    enddo
    do i=comp+1,this%N_mode
        Call lat%set_order_point(this%order(i),mode_base)
        tmp_internal(:,i)=mode_base(this%dat(i)%col)
    enddo
    vec=product(tmp_internal,dim=2)
    nullify(mode_base)
end subroutine

subroutine mode_reduce_comp_disc(this,ind_in,vec_in,comp,ind_out,vec_out)
    class(F_mode_rankN_sparse_col),intent(in)   :: this
    integer,intent(in)                          :: ind_in(:)
    real(8),intent(in)                          :: vec_in(:)
    integer,intent(in)                          :: comp 
    integer,intent(in)                          :: ind_out(:)
    real(8),intent(inout)                       :: vec_out(:)

    integer     :: i, j, i_out
    logical     :: mask(size(vec_in))

    do i=1,size(ind_out)
        i_out=ind_out(i)
        mask=.false.
        do j=1,this%ratio(comp)
            mask=mask.or.ind_in==this%dat(comp)%reverse(j,i_out)
        enddo
        vec_out(i)=sum(vec_in,mask=mask)
    enddo
end subroutine

subroutine mode_reduce_comp(this,lat,vec_in,comp,vec_out)
    use, intrinsic :: iso_fortran_env, only : error_unit
    class(F_mode_rankN_sparse_col),intent(in)   :: this
    real(8),intent(in)                          :: vec_in(:)
    type(lattice),intent(in)                    :: lat       !lattice type which knows about all states
    integer,intent(in)                          :: comp 
    real(8),intent(out)                         :: vec_out(lat%dim_modes(this%order(comp))*lat%Ncell)

    integer     ::  i

    vec_out=0.0d0
    do i=1,this%dat(comp)%nnz
        vec_out(this%dat(comp)%col(i))=vec_out(this%dat(comp)%col(i))+vec_in(i)
    enddo

    !!alternative way
    !do i=1,this%dat(comp)%dim_mat(2)
    !    vec_out(i)=sum(vec_in(this%dat(comp)%reverse(:,i)))
    !enddo
end subroutine

subroutine get_mode(this,lat,mode,tmp)
    class(F_mode_rankN_sparse_col),intent(in)   :: this
    type(lattice),intent(in)                    :: lat       !lattice type which knows about all states
    real(8),intent(out),pointer                 :: mode(:)   !pointer to required mode
    real(8),allocatable,target,intent(inout)    :: tmp(:)

    real(8)         :: tmp_internal(this%mode_size,size(this%dat))
    real(8),pointer :: mode_base(:)

    integer         :: i

    allocate(tmp(this%mode_size),source=0.0d0)
    mode=>tmp

    tmp_internal=1.d0
    do i=1,size(this%dat)
        Call lat%set_order_point(this%order(i),mode_base)
        tmp_internal(:,i)=mode_base(this%dat(i)%col)
    enddo
    mode=product(tmp_internal,dim=2)
    nullify(mode_base)
end subroutine

function is_same(this,comp)result(same)
    class(F_mode_rankN_sparse_col),intent(in)  :: this
    class(F_mode),intent(in)                   :: comp
    logical                                    :: same

    ERROR STOP "IMPLEMENT"
end function

subroutine destroy(this)
    class(F_mode_rankN_sparse_col),intent(inout) ::  this
    integer ::  i
    this%mode_size=-1
    this%order_occ=0
    do i=1,size(this%dat)
        deallocate(this%dat(i)%col)
        deallocate(this%dat(i)%reverse)
    enddo
    deallocate(this%dat)
    deallocate(this%order)
end subroutine

subroutine copy(this,F_out)
    class(F_mode_rankN_sparse_col),intent(in)    :: this
    class(F_mode),allocatable,intent(inout)     :: F_out

    integer ::  i

    Call this%copy_base(F_out)
    select type(F_out)
    class is(F_mode_rankN_sparse_col)
        F_out%mode_size=this%mode_size
        F_out%ratio=this%ratio
        allocate(F_out%dat(size(this%dat)))
        do i=1,size(this%dat)
            F_out%dat(i)%dim_mat=this%dat(i)%dim_mat
            F_out%dat(i)%nnz=this%dat(i)%nnz
            allocate(F_out%dat(i)%col,source=this%dat(i)%col)
            allocate(F_out%dat(i)%reverse,source=this%dat(i)%reverse)
        enddo
    class default
        ERROR STOP "FAILED TO COPY F_mode_rankN_sparse_col mode to F_out"
    end select
end subroutine

subroutine init_order(this,lat,abbrev_in,mat)
    use m_derived_types, only: op_abbrev_to_int
    class(F_mode_rankN_sparse_col),intent(inout) :: this
    type(lattice),intent(in)                :: lat       !lattice type which knows about all states
    character(len=*), intent(in)            :: abbrev_in !considered order abbreviations
    type(coo_mat),intent(inout)             :: mat(:)    !input matrices, destroyed when returned
    integer     :: order(len(abbrev_in))
    integer     :: i     

    order=op_abbrev_to_int(abbrev_in)
    Call this%init_base(order,mat(1)%dim_mat(1))
    if(size(mat)/=len(abbrev_in)) ERROR STOP "Matrix size has the be the same as the length of abbrev_in"
    allocate(this%dat(size(mat)))
    allocate(this%ratio(size(mat)),source=0)
    do i=1,size(mat)
        this%dat(i)%nnz     =mat(i)%nnz
        this%dat(i)%dim_mat =mat(i)%dim_mat
        this%ratio(i)=mat(i)%dim_mat(1)/mat(i)%dim_mat(2)
        Call move_alloc(mat(i)%col,this%dat(i)%col)
        Call mat(i)%destroy()
        Call set_reverse(this%dat(i))
    enddo
    this%mode_size=this%dat(1)%dim_mat(1)
end subroutine

subroutine set_reverse(mat)
    type(col_mat),intent(inout)     ::  mat
    integer :: i,j
    integer :: loc(mat%dim_mat(2))
    integer :: ratio

    ratio=mat%dim_mat(1)/mat%dim_mat(2)
    allocate(mat%reverse(ratio,mat%dim_mat(2)),source=0)
    loc=0
    do i=1,mat%nnz
        j=mat%col(i)
        loc(j)=loc(j)+1
        mat%reverse(loc(j),j)=i
    enddo
    if(minval(loc)/=ratio.or.maxval(loc)/=ratio) ERROR STOP "reverse col matrix failed, too compilcated state construction?"    !assumes every vector component appears ratio-times in the final state
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!            MPI ROUTINES           !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine bcast(this,comm)
    use mpi_basic                
    class(F_mode_rankN_sparse_col),intent(inout) ::  this        !this might fail if the server threads non-allocated class(F_mode), TAKE CARE OF THIS IN HAM_BASE
    type(mpi_type),intent(in)               ::  comm
#ifdef CPP_MPI
    integer     :: ierr, i
  
    Call this%bcast_base(comm)
    if(.not.allocated(this%ratio))then
        allocate(this%ratio(this%N_mode))
        allocate(this%dat(this%N_mode))
    endif
    Call MPI_Bcast(this%ratio,this%N_mode, MPI_INTEGER, comm%mas, comm%com,ierr)

    !bcast mat
    do i=1,this%N_mode
        Call this%dat(i)%bcast(comm)
    enddo
    Call MPI_Bcast(this%ratio, size(this%ratio), MPI_INTEGER, comm%mas, comm%com,ierr)
#else
    continue
#endif
end subroutine 

subroutine col_mat_bcast(this,comm)
    use mpi_util
    class(col_mat),intent(inout):: this
    type(mpi_type),intent(in)   :: comm

#ifdef CPP_MPI
    integer     :: ierr
    integer     :: shp(6)

    if(comm%ismas) shp=[this%dim_mat(1), this%dim_mat(2), this%nnz, size(this%col),size(this%reverse,1),size(this%reverse,2)]
    Call MPI_Bcast(shp, 6, MPI_INTEGER, comm%mas, comm%com,ierr)
    if(.not.comm%ismas)then
        this%dim_mat=shp(1:2)
        this%nnz=shp(3)
        if(.not.allocated(this%col)    ) allocate(this%col(shp(4)))
        if(.not.allocated(this%reverse)) allocate(this%reverse(shp(5),shp(6)))
    endif
    Call mpi_BCast(this%col,     size(this%col),     MPI_INTEGER, comm%mas, comm%com, ierr)
    Call mpi_BCast(this%reverse, size(this%reverse), MPI_INTEGER, comm%mas, comm%com, ierr)
#else
    continue
#endif
end subroutine

subroutine col_mat_send(this,ithread,tag,com)
    use mpi_util
    class(col_mat),intent(in)   :: this
    integer,intent(in)          :: ithread
    integer,intent(in)          :: tag
    integer,intent(in)          :: com

#ifdef CPP_MPI
    integer     :: ierr
    integer     :: shp(6)

    shp=[this%dim_mat(1), this%dim_mat(2), this%nnz, size(this%col),size(this%reverse,1),size(this%reverse,2)]
    Call mpi_send(shp,          6,                  MPI_INTEGER, ithread, tag, com, ierr)
    Call mpi_send(this%col,     size(this%col),     MPI_INTEGER, ithread, tag, com, ierr)
    Call mpi_send(this%reverse, size(this%reverse), MPI_INTEGER, ithread, tag, com, ierr)
#else
    continue
#endif
end subroutine

subroutine col_mat_recv(this,ithread,tag,com)
    use mpi_util
    class(col_mat),intent(inout):: this
    integer,intent(in)          :: ithread
    integer,intent(in)          :: tag
    integer,intent(in)          :: com

#ifdef CPP_MPI
    integer     :: ierr
    integer     :: shp(6)
    integer     :: stat(MPI_STATUS_SIZE) 

    Call mpi_recv(shp, 6, MPI_INTEGER, ithread, tag, com, stat, ierr)
    this%dim_mat=shp(1:2)
    this%nnz=shp(3)
    allocate(this%col(shp(4)),this%reverse(shp(5),shp(6)))
    Call mpi_recv(this%col,     size(this%col),     MPI_INTEGER, ithread, tag, com, stat, ierr)
    Call mpi_recv(this%reverse, size(this%reverse), MPI_INTEGER, ithread, tag, com, stat, ierr)
#else
    continue
#endif
end subroutine

subroutine send(this,ithread,tag,com)
    use mpi_basic                
    class(F_mode_rankN_sparse_col),intent(in)    :: this
    integer,intent(in)          :: ithread
    integer,intent(in)          :: tag
    integer,intent(in)          :: com

#ifdef CPP_MPI
    integer     :: ierr
    integer     :: N, i

    Call this%send_base(ithread,tag,com)
    N=size(this%dat)
    Call mpi_send(N, 1, MPI_INTEGER, ithread, tag, com, ierr)
    do i=1,N
        Call this%dat(i)%send(ithread,tag,com)
    enddo
    Call mpi_send(this%ratio, N, MPI_INTEGER, ithread, tag, com, ierr)
#else
    continue
#endif
end subroutine

subroutine recv(this,ithread,tag,com)
    use mpi_basic                
    class(F_mode_rankN_sparse_col),intent(inout) :: this
    integer,intent(in)          :: ithread
    integer,intent(in)          :: tag
    integer,intent(in)          :: com

#ifdef CPP_MPI
    integer     :: ierr
    integer     :: stat(MPI_STATUS_SIZE) 
    integer     :: N, i

    Call this%recv_base(ithread,tag,com)
    Call mpi_recv(N, 1, MPI_INTEGER, ithread, tag, com, stat, ierr)
    allocate(this%dat(N),this%ratio(N))
    do i=1,N
        Call this%dat(i)%recv(ithread,tag,com)
    enddo
    Call mpi_recv(this%ratio, N, MPI_INTEGER, ithread, tag, com, stat, ierr)
#else
    continue
#endif
end subroutine

end module
