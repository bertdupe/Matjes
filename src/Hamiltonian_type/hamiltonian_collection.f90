module m_hamiltonian_collection
use m_H_public, only: t_H, get_Htype_N
use m_H_type,only : len_desc
use m_derived_types, only: lattice
use mpi_basic
use, intrinsic  ::  ISO_FORTRAN_ENV, only: error_unit
private
public  ::  hamiltonian

type    ::  hamiltonian
    logical                     :: is_set=.false.
    class(t_H),allocatable      :: H(:)
    logical                     :: is_para(2)=.false. !signifies if any of the internal mpi-parallelizations have been initialized
    type(mpi_type)              :: com_global
    type(mpi_distv)             :: com_outer
    type(mpi_type)              :: com_inner
    integer                     :: NH_total=0  !global size of H
    integer                     :: NH_local=0  !local size of H (relevant with parallelization)

    character(len=len_desc),allocatable ::  desc_master(:)
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
    procedure   :: is_master
    procedure   :: bcast        !allows to bcast the Hamiltonian along one communicator before internal parallelization is done
    procedure   :: distribute   !distributes the different Hamiltonian entries to separate threads

    !small access routines
    procedure   :: size_H
    procedure   :: get_desc
end type

contains


subroutine distribute(this,com_in)
    use mpi_util
    use mpi_distrib_v
    !subroutine which distributes the H-Hamiltonians from the master of comm to have as few as possible H per thread
    !expects the hamiltonian of the comm-master to be initialized correctly
    class(hamiltonian),intent(inout)    :: this
    type(mpi_type),intent(in)           :: com_in
#ifdef CPP_MPI

    !integer         :: div, color
    integer         :: ierr

    integer         :: iH,ithread,i

    type(mpi_distv)             :: com_outer
    type(mpi_type)              :: com_inner

    class(t_H),allocatable      :: H_tmp(:)

    if(com_in%ismas)then
        if(.not.this%is_set) ERROR STOP "Cannot distribute hamiltonian that as not been initialized"
    endif
    Call bcast(this%NH_total,com_in)
    this%NH_local=this%NH_total
    if(com_in%Np==1) return
    !decide how to parallelize Hamiltonian
    Call get_two_level_comm(com_in,this%NH_total,com_outer,com_inner)

    if(com_outer%Np>1.and.com_inner%ismas) this%is_para(1)=.true.
    Call reduce_lor(this%is_para(1),com_in)

    !distribute the Hamiltonian to the masters of the inner parallelization
    if(com_inner%ismas)then
        if(com_outer%ismas)then
            !save descriptions for later easier io-access
            allocate(this%desc_master(this%NH_total))
            this%desc_master=this%H(:)%desc

            !send Hamiltonians
            do ithread=2,com_outer%Np
                do i=1,com_outer%cnt(ithread)
                    iH=com_outer%displ(ithread)+i
                    Call this%H(iH)%send(ithread-1,i,com_outer%com)
                enddo
            enddo
            !keep own ham entry
            this%NH_local=com_outer%cnt(1)
            Call get_Htype_N(H_tmp,this%NH_local)
            do i=1,this%NH_local
                Call this%H(i)%copy(H_tmp(i))    !moving would be much smarter, but has to be implemented
            enddo
            do i=1,size(this%H)
                Call this%H(i)%destroy()
            enddo
            deallocate(this%H)
            Call move_alloc(H_tmp,this%H)
        else
            this%NH_local=com_outer%cnt(com_outer%id+1)
            Call get_Htype_N(this%H,this%NH_local)
            do i=1,this%NH_local
                Call this%H(i)%recv(0,i,com_outer%com)
            enddo
        endif
        this%is_set=.true.
    endif

    if(com_inner%Np>1)then
        this%is_para(2)=.true.
        Call bcast(this%NH_local,com_inner)
        if(.not.com_inner%ismas) Call get_Htype_N(this%H,this%NH_local)
        do i=1,this%NH_local
            Call this%H(i)%distribute(com_inner)
        enddo
        this%is_set=.true.
    endif

    this%com_outer=com_outer
    this%com_inner=com_inner
    this%com_global=com_in
#else
    continue
#endif

end subroutine


subroutine bcast(this,comm)
    !bcast assuming Hamiltonian has not been scattered yet 
    class(hamiltonian),intent(inout)    :: this
    type(mpi_type),intent(in)           :: comm
#ifdef CPP_MPI
    integer     ::  i,N

    if(comm%ismas)then
        if(any(this%is_para))then
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
    if(this%is_para(1))then
        if(this%com_outer%ismas)then
            desc=this%desc_master(i)
        else
            ERROR STOP "CAN ONLY USE GET_DESC FROM MASTER THREAD"
        endif
    else
        desc=this%H(i)%desc
    endif
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
    this%NH_local=this%NH_total
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
    this%NH_local=this%NH_total
    this%is_set=.true.
end subroutine


subroutine energy_distrib(this,lat,order,Edist)
    !gets the energy at the sites, so far very wastefull with the memory
    use mpi_util
    use mpi_distrib_v
    class(hamiltonian),intent(inout)    :: this
    type(lattice), intent(in)           :: lat
    integer,intent(in)                  :: order
    real(8),allocatable,intent(inout)   :: Edist(:,:)

    integer     :: i
    integer     :: mult(this%com_outer%Np)
    
    if(allocated(Edist))then
        if(size(Edist,1)/=lat%Ncell*lat%site_per_cell(order).and.size(Edist,2)==this%NH_total) deallocate(Edist)
    endif
    if(.not.allocated(Edist)) allocate(Edist(lat%Ncell*lat%site_per_cell(order),this%NH_total),source=0.0d0)
    do i=1,this%NH_local
        Call this%H(i)%energy_dist(lat,order,Edist(:,i))
    enddo
    if(this%is_para(2)) Call reduce_sum(Edist,this%com_inner)
    if(this%is_para(1))then
        mult=0
        mult(this%com_outer%id+1)=size(Edist,1)*this%NH_local
        Call reduce_sum(mult,this%com_outer)
        Call gatherv(Edist,mult,this%com_outer)
    endif
end subroutine

subroutine energy_resolved(this,lat,E)
    use mpi_distrib_v, only: gatherv
    use mpi_util
    !get contribution-resolved energies
    !only correct on outer master thread
    class(hamiltonian),intent(in)   :: this
    type (lattice),intent(in)       :: lat
    real(8),intent(out)             :: E(this%NH_total)

    integer     ::  i

    E=0.0d0
    do i=1,this%NH_local
        Call this%H(i)%eval_all(E(i),lat)
    enddo
    if(this%is_para(2)) Call reduce_sum(E,this%com_inner)
    if(this%is_para(1).and.this%com_inner%ismas) Call gatherv(E,this%com_outer)
end subroutine

function energy(this,lat)result(E)
    !returns the total energy of the Hamiltonian array
    class(hamiltonian),intent(in)   :: this
    type(lattice),intent(in)        :: lat
    real(8)                         :: E

    real(8)     ::  tmp_E(this%NH_total)
    
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

    if(any(this%is_para)) ERROR STOP "IMPLEMENT"
    E=0.0d0
    do i=1,this%NH_total
        Call this%H(i)%eval_single(tmp_E(i),i_m,dim_bnd,lat)
    enddo
    tmp_E=tmp_E*real(this%H(:)%mult_M_single,8)
    E=sum(tmp_E)
end function

subroutine get_eff_field(this,lat,B,Ham_type,tmp)
    !calculates the effective internal magnetic field acting on the magnetization for the dynamics
    use mpi_basic
    use mpi_util
    class(hamiltonian),intent(in)       :: this
    type (lattice),intent(in)           :: lat    !lattice containing current order-parameters 
    real(8),intent(out)                 :: B(:)
    integer,intent(in)                  :: Ham_type   !integer that decides with respect to which mode the Hamiltonians derivative shall be obtained [1,number_different_order_parameters]
    real(8),intent(out)                 :: tmp(size(B))

    integer     :: iH, ierr
    B=0.0d0
    do iH=1,this%NH_local
        Call this%H(iH)%deriv(Ham_type)%get(this%H(iH),lat,B,tmp)
    enddo

    if(any(this%is_para)) Call reduce_sum(B,this%com_global)
    B=-B    !field is negative derivative
end subroutine

function is_master(this)result(master)
    use mpi_basic
    class(hamiltonian),intent(in)       :: this
    logical                             :: master

    integer :: ierr

    if(this%is_para(1))then
        master=this%com_outer%ismas
    else
        if(this%is_para(2))then
            master=this%com_inner%ismas
        else
            master=.true.
        endif
    endif
end function

end module
