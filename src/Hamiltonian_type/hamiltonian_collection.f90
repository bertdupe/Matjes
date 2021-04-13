module m_hamiltonian_collection
use m_H_public, only: t_H, get_Htype_N
use m_derived_types, only: lattice
use, intrinsic  ::  ISO_FORTRAN_ENV, only: error_unit
private
public  ::  hamiltonian

type    ::  hamiltonian
    logical                     :: is_set=.false.
    class(t_H),allocatable      :: H(:)
    logical                     :: para(2)=.false. !signifies if any of the internal mpi-parallelizations have been initialized
    integer                     :: NH_total=0  !global size of H
contains
    procedure   ::  init_H_mv
    procedure   ::  init_H_cp

    !energy getting routines
    procedure   :: energy
    procedure   :: energy_resolved
    procedure   :: energy_single
    procedure   :: energy_distrib
    
    !derivative getting routines
    procedure   :: get_eff_field

    !parallelization routines
    procedure   :: bcast    !allows to bcast the Hamiltonian along one communicator before internal parallelization is done

    !small access routines
    procedure   :: size_H
    procedure   :: get_desc
end type

contains

subroutine bcast(this,comm)
    !bcast assuming Hamiltonian has not been scattered yet 
    use mpi_basic                
    class(hamiltonian),intent(inout)    :: this
    type(mpi_type),intent(in)           :: comm
#ifdef CPP_MPI
    integer     ::  i,N

    if(comm%ismas)then
        if(any(this%para))then
            write(error_unit,'(3/A)') "Cannot broadcast Hamiltonian, since it appearst the Hamiltonian already has been scattered"  !world master only contains a part of the full Hamiltonian
            Error STOP
        endif
    endif

    Call MPI_Bcast(this%NH_total,1, MPI_INTEGER, comm%mas, comm%com,i)
    if(.not.comm%ismas) Call get_Htype_N(this%H,N)
    do i=1,this%NH_total
        Call this%H(i)%bcast(comm)
    enddo
#else
    continue
#endif
end subroutine

function size_H(this) result(N)
    class(Hamiltonian),intent(in)   :: this
    integer ::  N

    N=this%NH_total
end function

function get_desc(this,i) result(desc)
    use m_H_type, only: len_desc
    class(Hamiltonian),intent(in)       :: this
    integer,intent(in)                  :: i
    character(len=len_desc)             :: desc
    if(i<1.or.i>this%NH_total)then
        write(error_unit,'(3/A)') "Cannot get desciption of Hamiltonian, as the index is not with the bounds of the H-array"
        write(error_unit,'(A,I6)')  "Wanted index: ", i
        write(error_unit,'(A,I6)')  "H-arr size  : ", this%NH_total
        ERROR STOP
    endif
    desc=this%H(i)%desc
end function

subroutine init_H_mv(this,Harr)
    !initializes the Hamiltonian by moving the H array (thus destroying Harr)
    class(hamiltonian),intent(inout)        :: this
    class(t_H),allocatable,intent(inout)    :: Harr(:)
   
    if(.not.allocated(Harr))then
        write(error_unit,'(3/A)') "Cannot initialize Hamiltonian, since Harr-input is not allocated"
        ERROR STOP
    endif
    Call move_alloc(Harr,this%H)
    this%NH_total=size(this%H)
    this%is_set=.true.
end subroutine

subroutine init_H_cp(this,Harr)
    use, intrinsic  ::  ISO_FORTRAN_ENV, only: error_unit
    !initializes the Hamiltonian by moving the H array (thus destroying Harr)
    class(hamiltonian),intent(inout)        :: this
    class(t_H),allocatable,intent(in)       :: Harr(:)
    integer     ::  i
   
    if(.not.allocated(Harr))then
        write(error_unit,'(3/A)') "Cannot initialize Hamiltonian, since Harr-input is not allocated"
        ERROR STOP
    endif
    allocate(this%H,mold=Harr)
    do i=1,size(Harr)
        Call Harr(i)%copy(this%H(i))
    enddo
    this%NH_total=size(this%H)
    this%is_set=.true.
end subroutine


subroutine energy_distrib(this,lat,order,Edist)
    class(hamiltonian),intent(inout)    :: this
    type(lattice), intent(in)           :: lat
    integer,intent(in)                  :: order
    real(8),allocatable,intent(inout)   :: Edist(:,:)

    integer     ::  i
    
    if(allocated(Edist))then
        if(size(Edist,1)/=lat%Ncell*lat%site_per_cell(order).and.size(Edist,2)==this%NH_total) deallocate(Edist)
    endif
    if(.not.allocated(Edist)) allocate(Edist(lat%Ncell*lat%site_per_cell(order),this%NH_total),source=0.0d0)
    do i=1,this%NH_total
        Call this%H(i)%energy_dist(lat,order,Edist(:,i))
    enddo
end subroutine

subroutine energy_resolved(this,lat,E)
    !get contribution-resolved energies
    class(hamiltonian),intent(in)   :: this
    type (lattice),intent(in)       :: lat
    real(8),intent(out)             :: E(this%NH_total)

    integer     ::  i

    E=0.0d0
    do i=1,this%NH_total
        Call this%H(i)%eval_all(E(i),lat)
    enddo
end subroutine

function energy(this,lat)result(E)
    !returns the total energy of the Hamiltonian array
    class(hamiltonian),intent(in)   ::  this
    class(lattice),intent(in)       ::  lat
    real(8)                         ::  E

    real(8)     ::  tmp_E(size(this%H))
    
    Call this%energy_resolved(lat,tmp_E)
    E=sum(tmp_E)
end function

function energy_single(this,i_m,dim_bnd,lat)result(E)
    use m_derived_types, only: number_different_order_parameters
    !returns the total energy caused by a single entry !needs some updating 
    class(hamiltonian),intent(in)   :: this
    integer,intent(in)              :: i_m
    type (lattice),intent(in)       :: lat
    integer, intent(in)             :: dim_bnd(2,number_different_order_parameters)  !probably obsolete
    real(8)                         :: E

    real(8)     ::  tmp_E(this%NH_total)
    integer     ::  i

    E=0.0d0
    do i=1,this%NH_total
        Call this%H(i)%eval_single(tmp_E(i),i_m,dim_bnd,lat)
    enddo
    tmp_E=tmp_E*real(this%H(:)%mult_M_single,8)
    E=sum(tmp_E)
end function

subroutine get_eff_field(this,lat,B,Ham_type)
    !calculates the effective internal magnetic field acting on the magnetization for the dynamics
    class(hamiltonian),intent(inout)    :: this
    type (lattice),intent(in)           :: lat    !lattice containing current order-parameters 
    real(8),intent(inout)               :: B(:)
    integer,intent(in)                  :: Ham_type   !integer that decides with respect to which mode the Hamiltonians derivative shall be obtained [1,number_different_order_parameters]

    integer     :: iH
    real(8)     :: tmp(size(B))

    B=0.d0
    do iH=1,this%NH_total
        Call this%H(iH)%deriv(Ham_type)%get(this%H(iH),lat,B,tmp)
    enddo
    B=-B    !field is negative derivative
end subroutine


end module
