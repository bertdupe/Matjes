module m_dos_util
use m_TB_types, only: parameters_TB_IO_DOS
implicit none
private
public dos_nc, dos_sc, dos_bnd_nc, dos_bnd_sc, dos_orb_nc, dos_orb_sc, dos_t
real(8),parameter       ::  dist_inc=8.0d0  !how many sigma away from my the energy entries are still considered

type,abstract   :: dos_t
    private
    integer             :: N_entry=0    !number of times dos has been added !obsolete?, manually norm anyways...
    real(8)             :: sigma        !smearing sigma
    real(8),allocatable :: Eval(:)      !energy values
    real(8),allocatable :: dos(:)       !summed dos
    real(8)             :: dE
    real(8)             :: E_ext(2)
contains
    procedure :: init => init_dos
    procedure :: print=> print_dos
    procedure :: mpi_sum=> mpi_sum_dos
end type

type,extends(dos_t) ::  dos_nc
contains
    procedure :: add  => add_dos_nc
end type

type,extends(dos_t) :: dos_sc
contains
    procedure :: add  => add_dos_sc
end type


type,abstract,extends(dos_t)    ::  dos_bnd
    private
    integer     ::  bnd(2)
contains
    procedure :: init_bnd => init_dos_bnd
end type

type,extends(dos_bnd) ::  dos_bnd_nc
contains
    procedure :: add  => add_dos_bnd_nc
end type

type,extends(dos_bnd) ::  dos_bnd_sc
contains
    procedure :: add  => add_dos_bnd_sc
end type


type,abstract,extends(dos_t)    ::  dos_orb
    private
    integer     :: orb
    integer     :: freq
contains
    procedure :: init_bnd => init_dos_orb
end type

type,extends(dos_orb) ::  dos_orb_nc
contains
    procedure :: add  => add_dos_orb_nc
end type

type,extends(dos_orb) ::  dos_orb_sc
contains
    procedure :: add  => add_dos_orb_sc
end type

contains

subroutine init_dos_bnd(this,io_dos,bnd)
    class(dos_bnd),intent(inout)            :: this
    type(parameters_TB_IO_DOS),intent(in)   :: io_dos
    integer,intent(in)                      :: bnd(2)

    Call this%init(io_dos)
    this%bnd=bnd
end subroutine

subroutine init_dos_orb(this,io_dos,orb,ndim)
    class(dos_orb),intent(inout)            :: this
    type(parameters_TB_IO_DOS),intent(in)   :: io_dos
    integer,intent(in)                      :: orb, ndim

    Call this%init(io_dos)
    this%orb=orb
    this%freq=ndim
end subroutine

subroutine init_dos(this,io_dos)
    class(dos_t),intent(inout)              :: this
    type(parameters_TB_IO_DOS),intent(in)   :: io_dos

    integer             ::  NE,iE

    Ne=int((io_dos%E_ext(2)-io_dos%E_ext(1))/io_dos%dE)+1
    allocate(this%dos(Ne),source=0.0d0)
    allocate(this%Eval(Ne))
    do iE=1,Ne
        this%Eval(iE)=(iE-1)*io_dos%dE+io_dos%E_ext(1)
    enddo
    this%N_entry=0
    this%sigma=io_dos%sigma
    this%dE=io_dos%dE
    this%E_ext=io_dos%E_ext
end subroutine

subroutine add_dos_nc(this,eigval)
    class(dos_nc),intent(inout) :: this
    real(8),intent(in)          :: eigval(:)

    real(8)                     :: dos_loc(size(this%dos))
    integer                     :: i
    integer                     :: bnd(2)

    if(.not.allocated(this%Eval)) STOP "Trying to add to dos, by type is not initialized"
    dos_loc=0.0d0
    do i=1,size(eigval)
        Call get_ibnd (eigval(i),this,bnd)
        if(bnd(1)<bnd(2)) Call add_gauss(eigval(i),1.0d0,this%Eval(bnd(1):bnd(2)),dos_loc(bnd(1):bnd(2)),this%sigma)
    enddo
    this%dos=this%dos+dos_loc/real(size(eigval))
    this%N_entry=this%N_entry+1
end subroutine

subroutine add_dos_bnd_nc(this,eigval,eigvec)
    !weighted with the projections on a set contiguous set of orbitals (this%bnd(1):this%bnd(2)), irrespective of supercell
    class(dos_bnd_nc),intent(inout) :: this
    real(8),intent(in)              :: eigval(:)
    complex(8),intent(in)           :: eigvec(:,:)

    real(8)                     :: dos_loc(size(this%dos))
    integer                     :: i,bnd(2)
    real(8)                     :: pref 

    if(.not.allocated(this%Eval)) STOP "Trying to add to dos, by type is not initialized"
    dos_loc=0.0d0
    do i=1,size(eigval)
        pref=real(dot_product(eigvec(this%bnd(1):this%bnd(2),i),eigvec(this%bnd(1):this%bnd(2),i)),8)
        Call get_ibnd (eigval(i),this,bnd)
        if(bnd(1)<bnd(2)) Call add_gauss(eigval(i),pref,this%Eval(bnd(1):bnd(2)),dos_loc(bnd(1):bnd(2)),this%sigma)
    enddo
    this%dos=this%dos+dos_loc/real(size(eigval))
    this%N_entry=this%N_entry+1
end subroutine

subroutine add_dos_orb_nc(this,eigval,eigvec)
    !weighted with the projections on one orbital in the whole supercell (this%orb::this%freq)
    class(dos_orb_nc),intent(inout) :: this
    real(8),intent(in)              :: eigval(:)
    complex(8),intent(in)           :: eigvec(:,:)

    real(8)                     :: dos_loc(size(this%dos))
    integer                     :: i,bnd(2)
    real(8)                     :: pref 

    if(.not.allocated(this%Eval)) STOP "Trying to add to dos, by type is not initialized"
    dos_loc=0.0d0
    do i=1,size(eigval)
        pref=real(dot_product(eigvec(this%orb::this%freq,i),eigvec(this%orb::this%freq,i)),8)
        Call get_ibnd (eigval(i),this,bnd)
        if(bnd(1)<bnd(2)) Call add_gauss(eigval(i),pref,this%Eval(bnd(1):bnd(2)),dos_loc(bnd(1):bnd(2)),this%sigma)
    enddo
    this%dos=this%dos+dos_loc/real(size(eigval))
    this%N_entry=this%N_entry+1
end subroutine


subroutine add_dos_sc(this,eigval,eigvec)
    class(dos_sc),intent(inout) :: this
    real(8),intent(in)          :: eigval(:)
    complex(8),intent(in)       :: eigvec(:,:)

    real(8)                     :: dos_loc(size(this%dos))
    integer                     :: i, i_start, n_state, dimH, bnd(2)
    real(8)                     :: pref

    if(.not.allocated(this%Eval)) STOP "Trying to add to dos, by type is not initialized"
    !only consider the states above the zero energy
    i_start=0
    do i=1,size(eigval)
        if(eigval(i)>0.0d0)then
            i_start=i 
            exit
        endif
    enddo
    if(i_start==0) STOP "No positive eigenvalues found ind add_dos_sc. For dos choose energies in positive branch."

    dimH=size(eigvec,1)
    dos_loc=0.0d0
    do i=i_start,size(eigval)
        !u-part of BdG
        pref=real(dot_product(eigvec(1:dimH/2,i),eigvec(1:dimH/2,i)),8)
        Call get_ibnd (eigval(i),this,bnd)
        Call add_gauss(eigval(i),pref,this%Eval(bnd(1):bnd(2)),dos_loc(bnd(1):bnd(2)),this%sigma)
        !v-part of BdG
        pref=real(dot_product(eigvec(dimH/2+1:dimH,i),eigvec(dimH/2+1:dimH,i)),8)
        Call get_ibnd (-eigval(i),this,bnd)
        Call add_gauss(-eigval(i),pref,this%Eval(bnd(1):bnd(2)),dos_loc(bnd(1):bnd(2)),this%sigma)
    enddo
    this%dos=this%dos+dos_loc/real(size(eigval)-i_start+1)
    this%N_entry=this%N_entry+1
end subroutine

subroutine add_dos_bnd_sc(this,eigval,eigvec)
    class(dos_bnd_sc),intent(inout) :: this
    real(8),intent(in)          :: eigval(:)
    complex(8),intent(in)       :: eigvec(:,:)

    real(8)                     :: dos_loc(size(this%dos))
    integer                     :: i, i_start, n_state, dimH, bnd(2)
    real(8)                     :: pref

    if(.not.allocated(this%Eval)) STOP "Trying to add to dos, by type is not initialized"
    !only consider the states above the zero energy
    i_start=0
    do i=1,size(eigval)
        if(eigval(i)>0.0d0)then
            i_start=i 
            exit
        endif
    enddo
    if(i_start==0) STOP "No positive eigenvalues found ind add_dos_sc. For dos choose energies in positive branch."

    dimH=size(eigvec,1)
    dos_loc=0.0d0
    do i=i_start,size(eigval)
        !u-part of BdG
        pref=real(dot_product(eigvec(this%bnd(1):this%bnd(2),i),eigvec(this%bnd(1):this%bnd(2),i)),8)
        Call get_ibnd (eigval(i),this,bnd)
        if(bnd(1)<bnd(2)) Call add_gauss(eigval(i),pref,this%Eval(bnd(1):bnd(2)),dos_loc(bnd(1):bnd(2)),this%sigma)
        !v-part of BdG
        pref=real(dot_product(eigvec(dimH/2+this%bnd(1):dimH/2+this%bnd(2),i),eigvec(dimH/2+this%bnd(1):dimH/2+this%bnd(2),i)),8)
        Call get_ibnd (-eigval(i),this,bnd)
        if(bnd(1)<bnd(2)) Call add_gauss(-eigval(i),pref,this%Eval(bnd(1):bnd(2)),dos_loc(bnd(1):bnd(2)),this%sigma)
    enddo
    this%dos=this%dos+dos_loc/real(size(eigval)-i_start+1)
    this%N_entry=this%N_entry+1
end subroutine

subroutine add_dos_orb_sc(this,eigval,eigvec)
    class(dos_orb_sc),intent(inout) :: this
    real(8),intent(in)          :: eigval(:)
    complex(8),intent(in)       :: eigvec(:,:)

    real(8)                     :: dos_loc(size(this%dos))
    integer                     :: i, i_start, n_state, dimH, bnd(2)
    real(8)                     :: pref

    if(.not.allocated(this%Eval)) STOP "Trying to add to dos, by type is not initialized"
    !only consider the states above the zero energy
    i_start=0
    do i=1,size(eigval)
        if(eigval(i)>0.0d0)then
            i_start=i 
            exit
        endif
    enddo
    if(i_start==0) STOP "No positive eigenvalues found ind add_dos_sc. For dos choose energies in positive branch."

    dimH=size(eigvec,1)
    dos_loc=0.0d0
    do i=i_start,size(eigval)
        !u-part of BdG
        pref=real(dot_product(eigvec(this%orb:dimH/2:this%freq,i),eigvec(this%orb:dimH/2:this%freq,i)),8)
        Call get_ibnd (eigval(i),this,bnd)
        if(bnd(1)<bnd(2)) Call add_gauss(eigval(i),pref,this%Eval(bnd(1):bnd(2)),dos_loc(bnd(1):bnd(2)),this%sigma)
        !v-part of BdG
        pref=real(dot_product(eigvec(this%orb+dimH/2::this%freq,i),eigvec(this%orb+dimH/2::this%freq,i)),8)
        Call get_ibnd (-eigval(i),this,bnd)
        if(bnd(1)<bnd(2)) Call add_gauss(-eigval(i),pref,this%Eval(bnd(1):bnd(2)),dos_loc(bnd(1):bnd(2)),this%sigma)
    enddo
    this%dos=this%dos+dos_loc/real(size(eigval)-i_start+1)
    this%N_entry=this%N_entry+1
end subroutine


subroutine mpi_sum_dos(this,comm) 
    use mpi_util 
    class(dos_t),intent(inout)      :: this
    type(mpi_type),intent(in)       :: comm     !mpi communicator
    integer     :: ierr
#ifdef CPP_MPI

    Call reduce_sum(this%N_entry, comm)
    Call reduce_sum(this%dos, comm)
#else
    continue
#endif

end subroutine


subroutine add_gauss(val,pref,E,dos,sigma)
    !add to dos:
    !gauss distribution from single energy point val into the spacing supplied by E with std. sigma and an additional prefactor pref 
    use m_constants, only : pi
    real(8),intent(in)     ::  val,pref,sigma
    real(8),intent(in)     ::  E(:)
    real(8),intent(inout)  ::  dos(:)
    real(8)                ::  tmp(size(dos))
    tmp=(val-E)**2
    tmp=-tmp*0.5d0/sigma/sigma
    tmp=exp(tmp)
    tmp=tmp/sqrt(2.0d0*pi)/sigma
    dos=dos+pref*tmp
end subroutine

subroutine print_dos(this,fname)
    use, intrinsic :: iso_fortran_env, only : error_unit
    class(dos_t),intent(inout)  :: this
    character(len=*),intent(in) :: fname
    real(8),allocatable         :: dos_loc(:)
    real(8)                     :: norm
    integer ::  io,i
    if(.not.allocated(this%Eval)) STOP "Trying to print dos, by type is not initialized"
    if(this%N_entry<1) STOP "Trying to print dos, no eigenvalue sets have been added"
    open(newunit=io,file=fname)
    dos_loc=this%dos
    norm=sum(dos_loc)
    if(norm==0.0d0) then
        write(error_unit,'(3A)') "sum over all dos entries for ",fname," is 0, check dos input"
    else
        !dos_loc=dos_loc/norm
        dos_loc=dos_loc/real(this%N_entry,8)
    endif
    do i=1,size(this%dos)
       write(io,'(2E16.8)') this%Eval(i),dos_loc(i)
    enddo
    close(io)
end subroutine

!pure subroutine get_ibnd(val,dos,bnd)
subroutine get_ibnd(val,dos,bnd)
    !get the minimal-maximal index dist_inc*sigma around the considered energy eigenvalue needed for adding up the gauss distributions
    class(dos_t),intent(in) ::  dos
    real(8),intent(in)      ::  val
    integer,intent(out)     ::  bnd(2)
    bnd(1)=max(int(((val-dos%sigma*dist_inc)-dos%E_ext(1))/dos%dE)+1,1)
    bnd(2)=min(int(((val+dos%sigma*dist_inc)-dos%E_ext(1))/dos%dE),size(dos%Eval))
end subroutine
end module
