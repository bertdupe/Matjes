module m_work_ham_single
!module for derived type which provides some arrays for temporary array when evaluating the Energy of Hamiltonians
private
public work_ham_single

type work_ham_single
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
contains

subroutine finalize(this)
    type(work_ham_single),intent(inout)    ::  this

    Call this%destroy()
end subroutine

subroutine set(this,sizes)
    class(work_ham_single),intent(inout)    :: this
    integer,intent(in)                      :: sizes(2)

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
    class(work_ham_single),intent(inout)    :: this
    class(work_ham_single),intent(in)       :: arr(:)

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
    class(work_ham_single),intent(inout)    ::  this

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
