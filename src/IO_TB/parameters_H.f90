module m_parameters_TB_IO_H
!module of type which contains the integration kmesh part of the tight-binding
use m_input_H_types
use m_types_tb_h_inp
use m_ham_arrange
use m_delta_onsite
use m_wannier_inp, only: wann_dat

implicit none
private
public parameters_TB_IO_H
type parameters_TB_IO_H
    !parameters directly for the Hamiltonian
    !new parameters
    type(TB_hopping),allocatable                :: hop_io(:)        !hopping parameters
    type(TB_delta),allocatable                  :: del_io(:)        !delta parameters
    type(TB_Jsd),allocatable                    :: Jsd(:)           !delta parameters
    type(TBio_delta_onsite_scf),allocatable     :: del_scf_io(:)    !self-consistent delta on-site parameters
    type(TBio_defect),allocatable               :: defect(:)        !defect onsite-terms in super-cell

    type(Htb_inp)       ::  hop
    type(Hdelta)        ::  del

    !wann_io xor (wann_io_up and wann_io_dn) should be specified
    type(wann_dat)      ::  wann_io     !full wannier input data
    type(wann_dat)      ::  wann_io_up  !wannier input for up-states
    type(wann_dat)      ::  wann_io_dn  !wannier input for down-states

    real(8)             ::  Efermi=0.0d0   !Fermi energy added to the tight-binding hamiltonian
    integer             ::  nspin=1         !number of spins (1 or 2) for each orbital
    integer             ::  ncell=-1        !overall number of cells
    integer             ::  norb=-1         !number of orbitals in cell
    integer             ::  nsc=1           !2 if doubling for BdG superconductivity
    integer             ::  dimH=-1         !final size of Hamiltonian including all modifications
    integer,allocatable ::  norb_at(:)      !number of orbitals at each atom
    integer,allocatable ::  norb_at_off(:)  !offset of orbitals at each atom

    !solving parameters
    integer             ::  i_diag=1  !different diagonalization methods
    logical             ::  sparse=.false.  !do calculation sparse
    logical             ::  rearrange=.false.  !rearrange Hamiltonian basis order to have same site c and c^+  next to each other
    real(8)             ::  Ebnd(2)=[-1.0d+99,1.0d+99]     !minimal and maximal energy values to consider in restricted eigensolver routines
    integer             ::  estNe=0                       !estimated number of eigenvalues in interval
    real(8)             ::  diag_acc=1d-12    ! accuracy of iterative eigenvalue solution (so far only fpm input)

    !selfconsistent delta parameters
    logical             ::  use_scf=.false.         !use selfconsistent delta
    logical             ::  scf_print=.false.       !print intermediate delta steps
    integer             ::  scf_loopmax=100         !maximal number of loop iterations converging delta
    real(8)             ::  scf_diffconv=1.0d-6     !convergence criterion for difference of delta sum
    real(8)             ::  scf_Ecut=-1.0d0         !energy cutoff for selfconsistent delta energy sum 
    integer             ::  scf_kgrid(3)=[10,10,1]  !kgrid for delta-selfconsistency cycle in reciprocal space
    character(:),allocatable    ::  fname_kmesh     !file name for kmesh integration grid input

contains
    procedure   :: read_file
    procedure   :: bcast => bcast_local
end type
contains

subroutine bcast_local(this,comm) 
    use mpi_basic
    use mpi_util
    class(parameters_TB_IO_H),intent(inout)  :: this 
    type(mpi_type),intent(in)                :: comm
#ifdef CPP_MPI
    integer     :: i,N

#define __NAME__ hop_io
    N=0
    if(allocated(this%__NAME__)) N=size(this%__NAME__)
    Call bcast(N,comm)
    if(N>0.and..not.allocated(this%__NAME__)) allocate(this%__NAME__(N))
    do i=1,N
        Call this%__NAME__(i)%bcast(comm)
    enddo
#undef __NAME__

#define __NAME__ del_io
    N=0
    if(allocated(this%__NAME__)) N=size(this%__NAME__)
    Call bcast(N,comm)
    if(N>0.and..not.allocated(this%__NAME__)) allocate(this%__NAME__(N))
    do i=1,N
        Call this%__NAME__(i)%bcast(comm)
    enddo
#undef __NAME__

#define __NAME__ Jsd
    N=0
    if(allocated(this%__NAME__)) N=size(this%__NAME__)
    Call bcast(N,comm)
    if(N>0.and..not.allocated(this%__NAME__)) allocate(this%__NAME__(N))
    do i=1,N
        Call this%__NAME__(i)%bcast(comm)
    enddo
#undef __NAME__
    
#define __NAME__ del_scf_io
    N=0
    if(allocated(this%__NAME__)) N=size(this%__NAME__)
    Call bcast(N,comm)
    if(N>0.and..not.allocated(this%__NAME__)) allocate(this%__NAME__(N))
    do i=1,N
        Call this%__NAME__(i)%bcast(comm)
    enddo
#undef __NAME__

#define __NAME__ defect
    N=0
    if(allocated(this%__NAME__)) N=size(this%__NAME__)
    Call bcast(N,comm)
    if(N>0.and..not.allocated(this%__NAME__)) allocate(this%__NAME__(N))
    do i=1,N
        Call this%__NAME__(i)%bcast(comm)
    enddo
#undef __NAME__

    Call this%hop%bcast(comm)
    Call this%del%bcast(comm)

    Call this%wann_io%bcast(comm)
    Call this%wann_io_up%bcast(comm)
    Call this%wann_io_dn%bcast(comm)

    Call bcast      (this%Efermi,comm)
    Call bcast      (this%nspin,comm)
    Call bcast      (this%ncell,comm)
    Call bcast      (this%norb,comm)
    Call bcast      (this%nsc,comm)
    Call bcast      (this%dimH,comm)
    Call bcast_alloc(this%norb_at,comm)
    Call bcast_alloc(this%norb_at_off,comm)

    Call bcast      (this%i_diag,comm)
    Call bcast      (this%sparse,comm)
    Call bcast      (this%rearrange,comm)
    Call bcast      (this%Ebnd,comm)
    Call bcast      (this%estNe,comm)
    Call bcast      (this%diag_acc,comm)

    Call bcast      (this%use_scf,comm)
    Call bcast      (this%scf_print,comm)
    Call bcast      (this%scf_loopmax,comm)
    Call bcast      (this%scf_diffconv,comm)
    Call bcast      (this%scf_Ecut,comm)
    Call bcast      (this%scf_kgrid,comm)
    Call bcast_alloc(this%fname_kmesh,comm)

#else
    continue
#endif
end subroutine

subroutine read_file(this,io,fname)
    use m_io_read_util
    use m_io_utils
    use, intrinsic :: iso_fortran_env, only : output_unit, error_unit
    class(parameters_TB_IO_H),intent(inout)  :: this 
    integer,intent(in)                       :: io
    character(len=*), intent(in)             :: fname
    ! Internal variables
    integer :: N
    character(len=100)  :: str


    str=""
    call get_parameter(io,fname,'wann_ham',str)
    if(len_trim(str)/=0) Call this%wann_io%read_file(trim(adjustl(str)))

    if(.not.this%wann_io%is_set)then
        str=""
        call get_parameter(io,fname,'wann_up_ham',str)
        if(len_trim(str)/=0) Call this%wann_io_up%read_file(trim(adjustl(str)))

        str=""
        call get_parameter(io,fname,'wann_dn_ham',str)
        if(len_trim(str)/=0) Call this%wann_io_dn%read_file(trim(adjustl(str)))
    endif

    Call number_Hpar(io,fname,'TB_hopping',N)
    if(N>0)then
        allocate(this%hop_io(N))
        Call read_Hpar(io,fname,'TB_hopping',this%hop_io)
    elseif(.not.this%wann_io%is_set)then
        STOP "FAILED TO READ tight-binding hopping Hamiltonian (TB_hopping), but orbitals are set"
    endif

    Call number_Hpar(io,fname,'TB_delta',N)
    if(N>0)then
        allocate(this%del_io(N))
        Call read_Hpar(io,fname,'TB_delta',this%del_io)
    else
        write(output_unit,'(/A/)') "No tight-binding superconducting delta found"
    endif

    Call number_Hpar(io,fname,'TB_scfdelta',N)
    if(N>0)then
        allocate(this%del_scf_io(N))
        Call read_Hpar(io,fname,'TB_scfdelta',this%del_scf_io)
    else
        write(output_unit,'(/A/)') "No tight-binding self-consistent superconducting delta found"
    endif

    Call number_Hpar(io,fname,'TB_Jsd',N)
    if(N>0)then
        allocate(this%Jsd(N))
        Call read_Hpar(io,fname,'TB_Jsd',this%Jsd)
    else
        write(output_unit,'(/A/)') "No tight-binding Jsd-coupling found"
    endif

    Call number_Hpar(io,fname,'TB_defect',N)
    if(N>0)then
        allocate(this%defect(N))
        Call read_Hpar(io,fname,'TB_defect',this%defect)
    else
        write(output_unit,'(/A/)') "No tight-binding defect found"
    endif

    str=" "
    Call get_parameter(io,fname,'TB_scf_kmesh_file',str)
    if(len_trim(str)>1) this%fname_kmesh=trim(adjustl(str))

    call get_parameter(io,fname,'TB_Efermi',        this%Efermi)
    call get_parameter(io,fname,'TB_scf_print',     this%scf_print)
    call get_parameter(io,fname,'TB_scf_loopmax',   this%scf_loopmax)
    call get_parameter(io,fname,'TB_scf_diffconv',  this%scf_diffconv)
    call get_parameter(io,fname,'TB_scf_Ecut',      this%scf_Ecut)
    call get_parameter(io,fname,'TB_scf_kgrid',     this%scf_kgrid)
    call get_parameter(io,fname,'TB_sparse',        this%sparse)
    call get_parameter(io,fname,'TB_diag',          this%i_diag)
    call get_parameter(io,fname,'TB_diag_acc',      this%diag_acc)
    call get_parameter(io,fname,'TB_diag_Ebnd',     this%Ebnd)
    call get_parameter(io,fname,'TB_diag_Emin',     this%Ebnd(1))
    call get_parameter(io,fname,'TB_diag_Emax',     this%Ebnd(2))
    call get_parameter(io,fname,'TB_diag_estNe',    this%estNe)
    if(this%Ebnd(1)>=this%Ebnd(2))then
        write(error_unit,'(2/A/2(E16.8/))') "WARNING, tight binding minimal energy bound is smaller than maximal energy bound:", this%Ebnd
        STOP "Fix input"
    endif
end subroutine

subroutine number_Hpar(io,fname,var_name,Nnonzero)
    use m_io_read_util
    use, intrinsic :: iso_fortran_env, only : output_unit

    integer, intent(in)                         :: io
    character(len=*), intent(in)                :: fname,var_name
    integer,intent(out)                         :: Nnonzero

    ! internal variable
    integer :: Nentry
    integer :: i
    integer :: stat
    logical :: success
    class(TB_H_par),allocatable :: tmp
    character(len=100) :: str

    Nentry=0; Nnonzero=0
    Call set_pos_entry(io,fname,var_name,success)
    if(.not.success) return
    read(io,*) !read line with var_name, there is no additional information

    !We start to read the input
    !Find out how many entries there are 
    Call alloc_TB_H(tmp,var_name)
    do 
        read(io,'(a)',iostat=stat) str
        if (stat /= 0) STOP "UNEXPECTED END OF INPUT FILE"
        read(str,*,iostat=stat) tmp
        if (stat /= 0) exit
        Nentry=Nentry+1
        if(.not.tmp%is_zero()) Nnonzero=Nnonzero+1
    enddo
    Call write_info_number_found(Nentry,Nnonzero,var_name)
    success=.true.
    !rewind to read actual data
    do i=1,Nentry+1
        backspace(io)
    enddo
end subroutine

subroutine read_Hpar(io,fname,var_name,par)
    use m_io_read_util
    use, intrinsic :: iso_fortran_env, only : output_unit
    integer, intent(in)                         :: io
    character(len=*), intent(in)                :: fname,var_name
    class(TB_H_par), intent(inout)              :: par(:)

    ! internal variable
    integer :: ii
    integer :: stat
    character(len=100) :: str

    ii=1
    do while (ii<=size(par))
        read(io,'(a)',iostat=stat) str
        read(str,*,iostat=stat) par(ii)
        if(par(ii)%is_zero()) cycle
        write(output_unit,'(2A,I6,A)') var_name,' entry no.',ii,':'
        Call par(ii)%print_std()
        ii=ii+1
    enddo 
    Call check_further_entry(io,fname,var_name)
end subroutine

end module
