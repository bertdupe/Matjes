module m_mode_construction_rankN_sparse_col
use m_mode_construction
use m_type_lattice, only: dim_modes_inner,  lattice,number_different_order_parameters
use m_coo_mat
use m_work_ham_single, only:  work_mode, N_work
use, intrinsic  ::  ISO_FORTRAN_ENV, only: error_unit, output_unit
implicit none
private
public F_mode_rankN_sparse_col, col_mat

type col_mat
    integer :: dim_mat(2)=0
    integer :: nnz=0
    integer,allocatable :: col(:)
    integer,allocatable :: reverse(:,:)
contains
    procedure   :: is_same => col_mat_is_same
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
    procedure   :: get_mode_exc
    procedure   :: get_mode_exc_disc
    procedure   :: reduce_comp
    procedure   :: reduce_comp_add

    procedure   :: reduce_site_vec
    procedure   :: get_mode_single_size

    procedure   :: get_ind_site

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

    if(any(order==this%order))then
        order_ind=findloc(order==this%order,.true.,dim=1)
        dim_mode=dim_modes_inner(order)*this%ratio(order_ind)
    else
#ifdef CPP_DEBUG
        write(error_unit,'(A)') "Failed to get a single mode which is in the investigated order"
#endif
        dim_mode=0
    endif
end subroutine

subroutine get_ind_site(this,comp,site,size_out,ind)
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

subroutine get_mode_disc(this,lat,N,ind,vec)
    use, intrinsic :: iso_c_binding, only: C_PTR, C_loc
    use m_type_lattice, only: dim_modes_inner
    class(F_mode_rankN_sparse_col),intent(in)   :: this
    type(lattice),intent(in)                    :: lat
    integer,intent(in)                          :: N
    integer,intent(in)                          :: ind(N)
    real(8),intent(out)                         :: vec(N)

    real(8),pointer,contiguous :: mode_base(:)
    integer ::  i

    vec=1.0
    do i=1,size(this%dat)
        Call lat%set_order_point(this%order(i),mode_base)
        vec=vec*mode_base(this%dat(i)%col(ind))
    enddo
end subroutine

subroutine get_mode_exc_disc(this,lat,comp,N,ind,vec)
    use, intrinsic :: iso_c_binding, only: C_PTR, C_loc
    use m_type_lattice, only: dim_modes_inner
    class(F_mode_rankN_sparse_col),intent(in)   :: this
    type(lattice),intent(in)                    :: lat
    integer,intent(in)                          :: comp
    integer,intent(in)                          :: N
    integer,intent(in)                          :: ind(N)
    real(8),intent(out)                         :: vec(N)

    real(8),pointer,contiguous :: mode_base(:)
    integer ::  i

    vec=1.0
    do i=1,comp-1
        Call lat%set_order_point(this%order(i),mode_base)
        vec=vec*mode_base(this%dat(i)%col(ind))
    enddo
    do i=comp+1,this%N_mode
        Call lat%set_order_point(this%order(i),mode_base)
        vec=vec*mode_base(this%dat(i)%col(ind))
    enddo
end subroutine

subroutine get_mode_exc(this,lat,comp,work,vec)
    use, intrinsic :: iso_fortran_env, only : error_unit
    class(F_mode_rankN_sparse_col),intent(in)   :: this
    type(lattice),intent(in)                    :: lat       !lattice type which knows about all states
    integer,intent(in)                          :: comp
    type(work_mode),intent(inout)               :: work     !nothing to do for rank1
    real(8),intent(inout)                       :: vec(this%mode_size)

!    real(8)         :: tmp_modes(this%mode_size,this%N_mode-1)
    real(8),pointer,contiguous  :: tmp_modes(:,:)
    real(8),pointer,contiguous  :: mode_base(:)
    integer             :: i

    tmp_modes(1:this%mode_size,1:this%N_mode-1)=>work%real_arr(1+work%offset(1):this%mode_size*(this%N_mode-1)+work%offset(1))

    tmp_modes=1.d0
    do i=1,comp-1
        Call lat%set_order_point(this%order(i),mode_base)
        tmp_modes(:,i)=mode_base(this%dat(i)%col)
    enddo
    do i=comp+1,this%N_mode
        Call lat%set_order_point(this%order(i),mode_base)
        tmp_modes(:,i-1)=mode_base(this%dat(i)%col)
    enddo
    vec=product(tmp_modes,dim=2)
    nullify(mode_base)
end subroutine

subroutine reduce_comp(this,lat,vec_in,comp,vec_out)
    use, intrinsic :: iso_fortran_env, only : error_unit
    class(F_mode_rankN_sparse_col),intent(in)   :: this
    real(8),intent(in)                          :: vec_in(:)
    type(lattice),intent(in)                    :: lat       !lattice type which knows about all states
    integer,intent(in)                          :: comp 
    real(8),intent(out)                         :: vec_out(lat%dim_modes(this%order(comp))*lat%Ncell)

    integer     ::  i

    vec_out=0.0d0
    do i=1,this%dat(comp)%nnz
        vec_out(this%dat(comp)%col(i))=vec_in(i)
    enddo

    !!alternative way
    !do i=1,this%dat(comp)%dim_mat(2)
    !    vec_out(i)=sum(vec_in(this%dat(comp)%reverse(:,i)))
    !enddo
end subroutine

subroutine reduce_comp_add(this,lat,vec_in,comp,vec_out)
    use, intrinsic :: iso_fortran_env, only : error_unit
    class(F_mode_rankN_sparse_col),intent(in)   :: this
    real(8),intent(in)                          :: vec_in(:)
    type(lattice),intent(in)                    :: lat       !lattice type which knows about all states
    integer,intent(in)                          :: comp 
    real(8),intent(inout)                       :: vec_out(lat%dim_modes(this%order(comp))*lat%Ncell)

    integer     ::  i

    do i=1,this%dat(comp)%nnz
        vec_out(this%dat(comp)%col(i))=vec_out(this%dat(comp)%col(i))+vec_in(i)
    enddo

    !!alternative way
    !do i=1,this%dat(comp)%dim_mat(2)
    !    vec_out(i)=vec_out(i)+sum(vec_in(this%dat(comp)%reverse(:,i)))
    !enddo
end subroutine

subroutine get_mode(this,lat,mode,work,work_size)
    class(F_mode_rankN_sparse_col),intent(in)   :: this
    type(lattice),intent(in)                    :: lat       !lattice type which knows about all states
    real(8),intent(out),pointer,contiguous      :: mode(:)   !pointer to required mode
    type(work_mode),intent(inout)               :: work
    integer,intent(out)                         :: work_size(N_work)

    real(8),pointer,contiguous  :: tmp_modes(:,:)
    real(8),pointer,contiguous  :: mode_base(:)
    integer             :: i

    mode     (1:this%mode_size              )=>work%real_arr(1               +work%offset(1):this%mode_size                +work%offset(1))
    tmp_modes(1:this%mode_size,1:this%N_mode)=>work%real_arr(1+this%mode_size+work%offset(1):this%mode_size*(1+this%N_mode)+work%offset(1))

    tmp_modes=1.d0
    do i=1,size(this%dat)
        Call lat%set_order_point(this%order(i),mode_base)
        tmp_modes(:,i)=mode_base(this%dat(i)%col)
    enddo
    mode=product(tmp_modes,dim=2)
    nullify(mode_base,tmp_modes)
    work_size=0
    work_size(1)=this%mode_size
    work%offset=work%offset+work_size
end subroutine

function col_mat_is_same(this,comp)result(same)
    class(col_mat),intent(in)   :: this,comp
    logical                     :: same

    same=this%nnz==comp%nnz
    if(.not.same) return
    same=all(this%dim_mat==comp%dim_mat)
    if(.not.same) return
    same=allocated(this%col).and.allocated(comp%col)
    if(.not.same) return
    same=all(this%col==comp%col)
    if(.not.same) return
    same=allocated(this%reverse).and.allocated(comp%reverse)
    if(.not.same) return
    same=all(this%reverse==comp%reverse)
end function

function is_same(this,comp)result(same)
    class(F_mode_rankN_sparse_col),intent(in)  :: this
    class(F_mode                 ),intent(in)  :: comp
    logical                                    :: same

    integer     ::  i

    same=.false.
    select type (comp)
    class is(F_mode_rankN_sparse_col)
        same=all(this%ratio==comp%ratio)
        if(.not.same) return
        do i=1,size(this%dat)
            same=this%dat(i)%is_same(comp%dat(i))
            if(.not.same) return
        enddo
    end select
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
    class(F_mode_rankN_sparse_col),intent(in)   :: this
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
    class(F_mode_rankN_sparse_col),intent(inout)    :: this
    type(lattice),intent(in)                        :: lat       !lattice type which knows about all states
    character(len=*), intent(in)                    :: abbrev_in !considered order abbreviations
    type(coo_mat),intent(inout)                     :: mat(:)    !input matrices, destroyed when returned
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
