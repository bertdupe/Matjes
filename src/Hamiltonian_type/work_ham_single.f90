module m_work_ham_single
!module for derived type which provides some arrays for temporary array when evaluating the Energy of Hamiltonians
private
public work_ham_single, work_size_single, N_work_single, work_mode

integer,parameter   ::  N_work_single=2

type work_base
    !internal data arrays, DO NOT MODIFY UNLESS YOU KNOW WHAT YOU ARE DOING (pointers because targets in types does not work)
    real(8),pointer,contiguous, private :: real_dat(:)=>null()
    integer,pointer,contiguous, private :: int_dat(:)=>null()
    !data arrays to be used externally
    real(8),pointer,contiguous, public  :: real_arr(:)=>null()
    integer,pointer,contiguous, public  :: int_arr(:)=>null()
contains
    procedure   :: destroy
    procedure   :: set
    procedure   :: set_max
    final       :: finalize
end type

type,extends(work_base) ::  work_ham_single
end type

type,extends(work_base) ::  work_mode
end type


contains

subroutine work_size_single(dim_single,line_max,sizes)
    !routine which defines the size of the work array
    integer,intent(in)          :: dim_single   !corresponds to dim_single or dimT_single of the Hamiltonian implementation (dimension for inner work size array of set order)
    integer,intent(in)          :: line_max     !number of matrix entries on the col/row (corresponds to col_max/row_max)
    integer,intent(out)         :: sizes(N_work_single)

    sizes(1)=dim_single*(2*line_max+2)       !real array size
    sizes(2)=dim_single*(  line_max+2)+1     !integer array size
end subroutine


subroutine finalize(this)
    type(work_base),intent(inout)    ::  this

    Call this%destroy()
end subroutine

subroutine set(this,sizes)
    class(work_base),intent(inout)    :: this
    integer,intent(in)                :: sizes(2)

    Call this%destroy()
    if(sizes(1)>0)then
        allocate(this%real_dat(sizes(1)),source=0.0d0)
        this%real_arr=>this%real_dat
    endif
    if(sizes(2)>0)then
        allocate(this%int_dat (sizes(2)),source=0)
        this%int_arr=>this%int_dat
    endif
end subroutine

subroutine set_max(this,arr)
    class(work_base),intent(inout)    :: this
    class(work_base),intent(in)       :: arr(:)

    integer                                 :: sizes(2)
    integer                                 :: size_dat(size(arr))
    integer     ::  i
    
    do i=1,size(arr)
        size_dat(i)=size(arr(i)%real_dat)
    enddo
    sizes(1)=maxval(size_dat)

    do i=1,size(arr)
        size_dat(i)=size(arr(i)%int_dat)
    enddo
    sizes(2)=maxval(size_dat)

    Call this%set(sizes)
end subroutine

subroutine destroy(this)
    class(work_base),intent(inout)    ::  this

    if(associated(this%real_dat))then
        deallocate(this%real_dat)
        nullify(this%real_dat)
    endif
    if(associated(this%int_dat))then
        deallocate(this%int_dat)
        nullify(this%int_dat)
    endif
end subroutine
end module
