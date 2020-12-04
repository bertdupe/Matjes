module m_cell
use m_basic_types, only: atom
implicit none
private
public :: t_cell

type t_cell
    integer ::  n_attype
    type(atom), allocatable :: atomic(:)
    contains
    procedure :: ind_mag => get_ind_mag
    procedure :: num_mag => get_num_mag
    procedure :: get_magmom
    procedure :: get_mag_magmom
    procedure :: bcast
end type
contains

subroutine bcast(this,comm)
use mpi_basic                
    class(t_cell),intent(inout)    ::  this
    type(mpi_type),intent(in)       ::  comm

#ifdef CPP_MPI
    integer     :: ierr
    integer     :: N,i
    if(comm%ismas) N=size(this%atomic)
    Call MPI_Bcast(N, 1, MPI_INTEGER, comm%mas, comm%com,ierr)
    if(.not.comm%ismas) allocate(this%atomic(N))
    do i=1,N
        Call this%atomic(i)%bcast(comm)
    enddo   
#else
    continue
#endif
end subroutine

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

    Nat=this%num_mag()
    allocate(ind_Nat(Nat),source=0)
    j=0
    do i=1,size(this%atomic)
        if(this%atomic(i)%moment/=0.0d0)then
            j=j+1
            ind_Nat(j)=i
        endif
    enddo
end subroutine

function get_num_mag(this)result(nmag)
    class(t_cell),intent(in)    :: this
    integer                     :: nmag
    integer ::  i
    nmag=0
    do i=1,size(this%atomic)
        if(this%atomic(i)%moment/=0.0d0) nmag=nmag+1
    enddo
end function
    

end module
