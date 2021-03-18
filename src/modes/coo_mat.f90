module m_coo_mat
implicit none
private
public :: coo_mat
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
