module m_cell
use m_basic_types, only: atom
implicit none
private
public :: t_cell

type t_cell
     type(atom), allocatable :: atomic(:)
     contains
     procedure :: ind_mag => get_ind_mag
     procedure :: get_magmom
     procedure :: get_mag_magmom
end type
contains

subroutine get_magmom(this,magmom)
    !returns the magnetic moments of all atoms included in cell
    !(including the non-magnetic atoms)
    class(t_cell),intent(in)    ::  this
    real(8),allocatable,intent(out) ::  magmom(:)

    integer     ::  i,N

    N=size(this%atomic)
    allocate(magmom(N),source=0.0d0)
    do i=1,N
        magmom(i)=this%atomic(i)%moment
    enddo

end subroutine


subroutine get_mag_magmom(this,magmom)
    !returns the magnetic moments of all magnetic atoms included in cell
    class(t_cell),intent(in)    ::  this
    real(8),allocatable,intent(out) ::  magmom(:)

    integer,allocatable     ::  ind(:)

    integer     ::  i,N

    Call this%ind_mag(ind)
    N=size(ind)
    allocate(magmom(N),source=0.0d0)
    do i=1,N
        magmom(i)=this%atomic(ind(i))%moment
    enddo
    deallocate(ind)
end subroutine


subroutine get_ind_mag(this,ind_Nat)
    class(t_cell),intent(in)    ::  this
    integer                     ::  Nat
    integer                     ::  i,j
    integer,allocatable         ::  ind_Nat(:)

    Nat=0
    do i=1,size(this%atomic)
        if(this%atomic(i)%moment/=0.0d0) Nat=Nat+1
    enddo
    allocate(ind_Nat(Nat),source=0)
    j=0
    do i=1,size(this%atomic)
        if(this%atomic(i)%moment/=0.0d0)then
            j=j+1
            ind_Nat(j)=i
        endif
    enddo
end subroutine

end module
