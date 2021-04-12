module m_coo_mat
implicit none
private
public :: coo_mat, coo_full_unfold
type coo_mat
    integer :: dim_mat(2)=0
    integer :: nnz=0
    integer,allocatable :: row(:),col(:)
    real(8),allocatable :: val(:)
contains
    procedure   :: init
    procedure   :: copy
    procedure   :: destroy 
    procedure   :: mv
end type
contains

subroutine coo_full_unfold(rank,Ncell,dim_mode,mat)
    integer,intent(in)              :: rank
    integer,intent(in)              :: Ncell
    integer,intent(in)              :: dim_mode(rank)
    type(coo_mat),intent(out)       :: mat(rank)

    real(8),allocatable     :: val(:)
    integer,allocatable     :: row(:),col(:)
    integer     :: div
    integer     :: Nmode
    integer     :: dim_mode_full

    integer     :: i_site,i,i_mode,ind_site
    integer     :: ind,ii


    dim_mode_full=product(dim_mode)
    Nmode=dim_mode_full*Ncell

    allocate(val(Nmode),source=1.0d0)
    allocate(row(Nmode),col(Nmode))
    row=[(i,i=1,Nmode)]
    do i_mode=1,rank
        ii=0
        div=product(dim_mode(:i_mode-1))
        do i_site=1,Ncell
            ind_site=(i_site-1)*dim_mode(i_mode)
            do i=1,dim_mode_full
                ind=(i-1)/div
                ind=modulo(ind,dim_mode(i_mode))+1+ind_site
                ii=ii+1
                col(ii)=ind
            enddo
        enddo
        Call mat(i_mode)%init([Nmode,dim_mode(i_mode)],Nmode,row,col,val)
    enddo
    deallocate(val,row,col)

end subroutine


subroutine init(this,dim_mat,nnz,row,col,val)
    class(coo_mat),intent(inout)    :: this
    integer,intent(in)              :: dim_mat(2),nnz
    integer,intent(in)              :: row(nnz),col(nnz)
    real(8),intent(in)              :: val(nnz)

    if(allocated(this%row)) deallocate(this%row)
    if(allocated(this%col)) deallocate(this%col)
    if(allocated(this%val)) deallocate(this%val)
    this%dim_mat=dim_mat
    this%nnz=nnz
    allocate(this%row,source=row)
    allocate(this%col,source=col)
    allocate(this%val,source=val)
end subroutine

subroutine mv(this,mat_out)
    class(coo_mat),intent(inout)    :: this
    type(coo_mat),intent(inout)     :: mat_out

    if(allocated(mat_out%row)) deallocate(mat_out%row)
    if(allocated(mat_out%col)) deallocate(mat_out%col)
    if(allocated(mat_out%val)) deallocate(mat_out%val)
    Call move_alloc(this%row,mat_out%row)
    Call move_alloc(this%col,mat_out%col)
    Call move_alloc(this%val,mat_out%val)
    mat_out%dim_mat=this%dim_mat
    mat_out%nnz    =this%nnz
    Call this%destroy()
end subroutine

subroutine copy(this,mat_out)
    class(coo_mat),intent(in)       :: this
    type(coo_mat),intent(inout)     :: mat_out

    if(mat_out%nnz/=this%nnz)then
        if(allocated(mat_out%row)) deallocate(mat_out%row)
        if(allocated(mat_out%col)) deallocate(mat_out%col)
        if(allocated(mat_out%val)) deallocate(mat_out%val)
        allocate(mat_out%row(this%nnz))
        allocate(mat_out%col(this%nnz))
        allocate(mat_out%val(this%nnz))
    endif
    mat_out%nnz    =this%nnz
    mat_out%dim_mat=this%dim_mat
    mat_out%row    =this%row
    mat_out%col    =this%col
    mat_out%val    =this%val
end subroutine

subroutine destroy(this)
    class(coo_mat),intent(inout)    :: this

    if(allocated(this%row)) deallocate(this%row)
    if(allocated(this%col)) deallocate(this%col)
    if(allocated(this%val)) deallocate(this%val)
    this%nnz=0
    this%dim_mat=0
end subroutine
end module
