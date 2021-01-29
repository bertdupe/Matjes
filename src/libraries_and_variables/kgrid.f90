module m_kgrid
use m_type_lattice, only : lattice
private
public k_grid_t

type k_grid_t
    integer             :: k_offset(3)
    integer             :: kgrid(3)
    real(8)             :: kdiff(3,3)
    integer             :: Nk=0

contains
    procedure :: set => k_grid_set
    procedure :: get_k  => k_grid_get_k
    procedure :: get_Nk  => k_grid_get_Nk
end type

contains
function k_grid_get_Nk(this)result(Nk)
    class(k_grid_t),intent(in)  :: this
    integer                     :: Nk
    Nk=this%Nk
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
