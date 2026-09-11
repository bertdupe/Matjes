module m_cell
use m_basic_types, only: atom
use, intrinsic  ::  ISO_FORTRAN_ENV, only: error_unit, output_unit
implicit none
private
public :: t_cell

type t_cell
    integer ::  n_attype=0
    type(atom), allocatable :: atomic(:)
    contains
    procedure :: ind_mag_all
    procedure :: ind_mag
    procedure :: ind_ph
    procedure :: ind_attype
    procedure :: get_magmom
    procedure :: get_mag_magmom
    procedure :: ind_Z_all
    procedure :: get_Z_phonon
    procedure :: ind_M_all
    procedure :: get_M_phonon
    procedure :: index_all_magnetic
    procedure :: type_all_magnetic
    procedure :: pos_all_magnetic
    procedure :: bcast
    procedure :: copy
end type
contains

subroutine pos_all_magnetic(this,pos)
    !subroutine to return the position of all the magnetic atoms in the cell
    class(t_cell),intent(in)            :: this
    real(8),allocatable,intent(inout)   :: pos(:)

    integer ::  i,nmag,j
    real(8) :: min_mag=1.0d-8


    if(all(this%atomic(:)%moment.lt.min_mag))then
        write(output_unit,'(a)') 'Failed to get atomic indices of magnetic atoms'
        STOP "MISTAKE SETTING UP HAMILTONIAN OR CELL?"
    endif
    if(allocated(pos)) deallocate(pos)
    nmag=count(this%atomic(:)%moment.gt.min_mag)
    allocate(pos(3*nmag),source=0.0d0)
    j=1
    do i=1,size(this%atomic)
       if(this%atomic(i)%moment.gt.min_mag) then
       pos(j:j+2)=this%atomic(i)%position
       j=j+3
       endif
    enddo
end subroutine

subroutine type_all_magnetic(this,ind)
    !subroutine to return the types of all the magnetic atoms in the cell
    class(t_cell),intent(in)            :: this
    integer,allocatable,intent(inout)   :: ind(:)

    integer ::  i
    real(8) :: min_mag=1.0d-8


    if(all(this%atomic(:)%moment.lt.min_mag))then
        write(output_unit,'(a)') 'Failed to get atomic indices of magnetic atoms'
        STOP "MISTAKE SETTING UP HAMILTONIAN OR CELL?"
    endif
    if(allocated(ind)) deallocate(ind)
    ind=pack([(this%atomic(i)%type_id,i=1,size(this%atomic))],this%atomic(:)%moment.gt.min_mag)
end subroutine

subroutine index_all_magnetic(this,ind)
    !subroutine to return the indices of the magnetic atoms in the cell
    class(t_cell),intent(in)            :: this
    integer,allocatable,intent(inout)   :: ind(:)

    integer ::  i
    real(8) :: min_mag=1.0d-8


    if(all(this%atomic(:)%moment.lt.min_mag))then
        write(output_unit,'(a)') 'Failed to get atomic indices of magnetic atoms'
        STOP "MISTAKE SETTING UP HAMILTONIAN OR CELL?"
    endif
    if(allocated(ind)) deallocate(ind)
    ind=pack([(i,i=1,size(this%atomic))],this%atomic(:)%moment.gt.min_mag)
end subroutine

subroutine ind_attype(this,id,ind)
    !subroutine to return the indices of the atoms in the cell which have the same type_id as the id-input
    class(t_cell),intent(in)            :: this
    integer,intent(in)                  :: id
    integer,allocatable,intent(inout)   :: ind(:)

    integer ::  i

    
    if(id<1.or.id>maxval(this%atomic(:)%type_id))then
        write(output_unit,'(a,I10)') 'Failed to get atomic indices of atom type',id
        STOP "MISTAKE SETTING UP HAMILTONIAN OR CELL?"
    endif
    if(allocated(ind)) deallocate(ind)
    ind=pack([(i,i=1,size(this%atomic))],this%atomic(:)%type_id==id)
end subroutine
    
subroutine copy(this,cell_out)
    class(t_cell),intent(in)        ::  this
    class(t_cell),intent(inout)     ::  cell_out

    cell_out%n_attype=this%n_attype
    allocate(cell_out%atomic,source=this%atomic)
end subroutine

subroutine bcast(this,comm)
    use mpi_basic                
    class(t_cell),intent(inout)    ::  this
    type(mpi_type),intent(in)       ::  comm

#ifdef CPP_MPI
    integer     :: ierr
    integer     :: N,i
    if(comm%ismas) N=size(this%atomic)
    Call MPI_Bcast(N, 1, MPI_INTEGER, comm%mas, comm%com,ierr)
    Call MPI_Bcast(this%n_attype, 1, MPI_INTEGER, comm%mas, comm%com,ierr)
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

    Call this%ind_mag_all(ind)
    N=size(ind)
    allocate(magmom(N),source=0.0d0)
    do i=1,N
        magmom(i)=this%atomic(ind(i))%moment
    enddo
    deallocate(ind)
end subroutine

subroutine get_Z_phonon(this,phonon)
    !returns the magnetic moments of all magnetic atoms included in cell
    class(t_cell),intent(in)    ::  this
    real(8),allocatable,intent(out) ::  phonon(:)

    integer,allocatable     ::  ind(:)

    integer     ::  i,N

    Call this%ind_Z_all(ind)
    N=size(ind)
    allocate(phonon(N),source=0.0d0)
    do i=1,N
        phonon(i)=this%atomic(ind(i))%charge
    enddo
    deallocate(ind)
end subroutine

subroutine get_m_phonon(this,phonon)
    !returns the masses of all atoms included in cell
    class(t_cell),intent(in)    ::  this
    real(8),allocatable,intent(out) ::  phonon(:)

    integer,allocatable     ::  ind(:)

    integer     ::  i,N

    Call this%ind_M_all(ind)
    N=size(ind)
    allocate(phonon(N),source=0.0d0)
    do i=1,N
        phonon(i)=this%atomic(ind(i))%mass
    enddo
    deallocate(ind)
end subroutine

function ind_ph(this,ind_at)
    class(t_cell),intent(in)    ::  this
    integer,intent(in)          ::  ind_at
    integer                     ::  ind_ph
    
    if(ind_at<1.or.ind_at>size(this%atomic)) ERROR STOP "Trying to get phonon index of atom outside of set atom range"
    if(.not.this%atomic(ind_at)%use_ph) ERROR STOP "Trying to get phonon index of atom with disabled phonons"
    ind_ph=count(this%atomic(:ind_at)%use_ph)
end function


function ind_mag(this,ind_at)
    !returns the magnetic index of a atom in the cell
    ! atom index -> xth magnetic atom
    class(t_cell),intent(in)    ::  this
    integer,intent(in)          ::  ind_at
    integer                     ::  ind_mag
    
    if(ind_at<1.or.ind_at>size(this%atomic)) ERROR STOP "Trying to get magnetic moment index of atom outside of set atom range"
    if(this%atomic(ind_at)%moment==0.0d0) ERROR STOP "Trying to get magnetic moment index of non-magnetic atom"
    ind_mag=count(this%atomic(:ind_at)%moment/=0.0d0)
end function

subroutine ind_mag_all(this,ind_Nat)
    !return the atomic indices of all magnetic atoms in the cell
    class(t_cell),intent(in)    ::  this
    integer,allocatable         ::  ind_Nat(:)
    integer     ::  i

    ind_Nat=pack([(i,i=1,size(this%atomic))],this%atomic(:)%moment/=0.0d0)
end subroutine

subroutine ind_Z_all(this,ind_Nat)
    class(t_cell),intent(in)    ::  this
    integer,allocatable         ::  ind_Nat(:)
    integer     ::  i

    ind_Nat=pack([(i,i=1,size(this%atomic))],this%atomic(:)%charge/=0.0d0)
end subroutine

subroutine ind_M_all(this,ind_Nat)
    class(t_cell),intent(in)    ::  this
    integer,allocatable         ::  ind_Nat(:)
    integer     ::  i

    ind_Nat=pack([(i,i=1,size(this%atomic))],this%atomic(:)%mass/=0.0d0)

end subroutine

end module
