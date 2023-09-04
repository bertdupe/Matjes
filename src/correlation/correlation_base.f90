module m_correlation_base
use m_io_files_utils
use m_convert
use, intrinsic :: iso_fortran_env, only : error_unit,output_unit
implicit none

private
public :: corre_base

type, abstract :: corre_base
    real(8),dimension(:),allocatable :: data                ! data to be correlated
    real(8),dimension(:),allocatable :: correlation         ! correlations on the data
    logical :: is_set=.false.                               ! do you calculate the correlations
    logical :: is_print=.false.                             ! print data
    integer,dimension(:),allocatable :: shape_data
    character(len=100)               :: name=""             ! name of the correlation. for example M(t,t') or M(t)*g(t-t')*M(t')


    ! defered type
    procedure(int_get_cor),pointer,nopass,public    :: get_correlation => null()

    contains

    ! unchanged stuff
    procedure, NON_OVERRIDABLE     :: init_base
    procedure, NON_OVERRIDABLE     :: copy
    procedure, NON_OVERRIDABLE     :: print_base
    procedure, NON_OVERRIDABLE     :: destroy_base
    procedure, NON_OVERRIDABLE     :: overwrite_base

end type

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
abstract interface

   pure subroutine int_get_cor(data,t,tprime,chunk,r)
      real(8),intent(in)    :: data(:)
      integer,intent(in)    :: t,tprime,chunk(:)
      real(8),intent(out)   :: r(:)
   end subroutine


end interface

contains

subroutine init_base(this,shape_data)
    class(corre_base),intent(inout) :: this
    integer,intent(in)              :: shape_data(:)

    integer :: N

    if (this%is_set) stop "correlations are already set"

    N=product(shape_data)
    allocate(this%data(N),source=0.0d0)
    allocate(this%correlation(N),source=0.0d0)

    N=size(shape_data)
    allocate(this%shape_data(N),source=shape_data)

    this%is_set=.true.

end subroutine

subroutine copy(this,out)
    class(corre_base),intent(in)  :: this
    class(corre_base),intent(out) :: out

    integer :: N

    if (this%is_set) stop "correlations are already set and cannot be copied"

    N=size(this%data)
    allocate(out%data(N),source=this%data)
    allocate(out%correlation(N),source=this%correlation)

    N=size(this%shape_data)
    allocate(out%shape_data(N),source=this%shape_data)

    out%name=this%name
    if (.not.associated(this%get_correlation)) stop "get_correlation is not associated and cannot be copied"
    out%get_correlation => this%get_correlation

    out%is_set=.true.

end subroutine


subroutine print_base(this,tag)
    class(corre_base),intent(in) :: this
    integer,intent(in)           :: tag

    integer :: i,io_out,offset
    character(len=100) :: fname,form

    if (.not.this%is_print) return

    fname=convert('corre_',tag)
    form=convert('(',this%shape_data(1),'(E20.12E3,2x))')
    offset=this%shape_data(1)

    io_out=open_file_write(trim(fname))

    do i=1,this%shape_data(2),offset
       write(io_out,form) this%correlation(i:i+offset)
    enddo

    call close_file(fname,io_out)

end subroutine

subroutine overwrite_base(this,x)
    class(corre_base),intent(inout) :: this
    real(8),intent(in)              :: x(:)

    if (size(x).ne.size(this%data)) stop "can not overwrite data - sizes do not match"
    this%data=x

end subroutine




subroutine destroy_base(this)
    class(corre_base),intent(inout) :: this

    if (.not.this%is_set) stop "correlations are not set"

    deallocate(this%data,this%correlation,this%shape_data)

    this%is_set=.false.

end subroutine

end module
