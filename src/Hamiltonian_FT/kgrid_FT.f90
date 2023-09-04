module m_kgrid_FT
use m_kgrid
implicit none

private
public :: k_grid_FT

type,extends(kmesh_t) :: k_grid_FT
    integer             :: k_offset(3)
    integer             :: kgrid(3)
    real(8)             :: kdiff(3,3)
    integer             :: Nk=0

contains
    procedure :: kmesh_read
    procedure :: get_kgrid
    procedure :: set            => k_grid_set
    procedure :: get_k          => k_grid_get_k
    procedure :: get_Nk         => k_grid_get_Nk
    procedure :: get_normalize  => k_grid_get_normalize
    procedure :: bcast          => k_grid_bcast
end type

contains

subroutine kmesh_read(this,prefix,io,fname)
    use m_convert
    use m_io_read_util
    use m_io_utils
    class(k_grid_FT),intent(inout)              :: this
    integer,intent(in)                          :: io
    character(len=*), intent(in)                :: fname,prefix

    this%kgrid=(/10,10,10/)
    this%k_offset=(/0.0d0,0.0d0,0.0d0/)


    call get_parameter(io,fname,prefix//'_kgrid',this%kgrid)

    this%Nk=product(this%kgrid)

end subroutine

function k_grid_get_Nk(this)result(Nk)
    class(k_grid_FT),intent(in)  :: this
    integer                     :: Nk
    Nk=this%Nk
end function

function get_kgrid(this)result(grid)
    class(k_grid_FT),intent(in)  :: this
    integer                     :: grid(3)
    grid=this%kgrid
end function

function k_grid_get_normalize(this)result(norm)
    class(k_grid_FT),intent(in)    :: this
    real(8)                       :: norm
    norm=1.0d0/real(this%Nk)
end function


subroutine k_grid_set(this,a_inv,kgrid)
    class(k_grid_FT),intent(out)     :: this
    real(8),intent(in)              :: a_inv(3,3)   !reciprocal space lattice vectors WITH 2*pi factor
    integer,intent(in)              :: kgrid(3)
    integer ::  i

    this%kgrid=kgrid
    this%Nk=product(this%kgrid)
    this%kdiff=a_inv/spread(real(this%kgrid),2,3)
    this%k_offset=[(product(this%kgrid(1:i)),i=0,2)]
end subroutine

function k_grid_get_k(this,i)result(k)
    use m_constants,    only : pi
    class(k_grid_FT),intent(in)      :: this
    integer,intent(in)              :: i
    real(8)                         :: k(3)
    integer                         :: kint(3)

    if(i<1.or.i>this%NK) ERROR STOP "wanted k-grid index not within kgrid-lattice"
    kint=modulo((i-1)/this%k_offset,this%kgrid)
    k=matmul(kint,this%kdiff)
end function

subroutine k_grid_bcast(this,comm)
    use mpi_util,only: bcast, mpi_type
    class(k_grid_FT),intent(inout)   :: this
    type(mpi_type),intent(in)       :: comm
#ifdef CPP_MPI
    integer     :: int_arr(7)

    int_arr(1)  =this%Nk
    int_arr(2:4)=this%kgrid
    int_arr(5:7)=this%k_offset
    Call bcast(int_arr,comm)
    Call bcast(this%kdiff,comm)
    this%Nk      =int_arr(1  )
    this%kgrid   =int_arr(2:4)
    this%k_offset=int_arr(5:7)
#endif
end subroutine



end module m_kgrid_FT
