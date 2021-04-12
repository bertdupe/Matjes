module m_kgrid
use m_type_lattice, only : lattice
private
public k_grid_t, kmesh_t

type,abstract   :: kmesh_t
contains
    procedure(int_get_Nk),deferred          :: get_nK
    procedure(int_get_k),deferred           :: get_K
    procedure(int_get_normalize),deferred   :: get_normalize
end type

type,extends(kmesh_t) :: k_grid_t
    integer             :: k_offset(3)
    integer             :: kgrid(3)
    real(8)             :: kdiff(3,3)
    integer             :: Nk=0

contains
    procedure :: set            => k_grid_set
    procedure :: get_k          => k_grid_get_k
    procedure :: get_Nk         => k_grid_get_Nk
    procedure :: get_normalize  => k_grid_get_normalize
end type

abstract interface
    function int_get_Nk(this)result(Nk)
        import kmesh_t
        class(kmesh_t),intent(in)  :: this
        integer                    :: Nk
    end function

    function int_get_k(this,i)result(k)
        import kmesh_t
        class(kmesh_t),intent(in)  :: this
        integer,intent(in)          :: i
        real(8)                     :: k(3)
    end function

    function int_get_normalize(this)result(norm)
        import kmesh_t
        class(kmesh_t),intent(in)    :: this
        real(8)                      :: norm
    end function

end interface


contains
function k_grid_get_Nk(this)result(Nk)
    class(k_grid_t),intent(in)  :: this
    integer                     :: Nk
    Nk=this%Nk
end function

function k_grid_get_normalize(this)result(norm)
    class(k_grid_t),intent(in)    :: this
    real(8)                       :: norm
    norm=1.0d0/real(this%Nk)
end function

    
subroutine k_grid_set(this,a_inv,kgrid)
    use m_constants,    only : pi
    class(k_grid_t),intent(out)     :: this
    real(8),intent(in)              :: a_inv(3,3)   !reciprocal space lattice vectors without 2*pi factor(so that the ofter used lat%a_sc_inv works directly)
    integer,intent(in)              :: kgrid(3)     
    integer ::  i
    
    this%kgrid=kgrid
    this%Nk=product(this%kgrid)
    this%kdiff=a_inv/spread(real(this%kgrid),2,3)*2.0d0*pi
    this%k_offset=[(product(this%kgrid(1:i)),i=0,2)]
end subroutine

function k_grid_get_k(this,i)result(k)
    use m_constants,    only : pi
    class(k_grid_t),intent(in)      :: this
    integer,intent(in)              :: i
    real(8)                         :: k(3)
    integer                         :: kint(3)

    if(i<1.or.i>this%NK) ERROR STOP "wanted k-grid index not within kgrid-lattice"
    kint=modulo((i-1)/this%k_offset,this%kgrid)
    k=matmul(kint,this%kdiff)
end function
end module 
