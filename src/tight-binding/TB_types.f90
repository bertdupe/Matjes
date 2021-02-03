module m_TB_types
use m_input_H_types
use m_types_tb_h_inp
use m_ham_arrange
use m_delta_onsite
use m_dos_io
use m_wannier_inp, only: wann_dat
implicit none
public
private upd_h_par, init_ham_init

type parameters_ham_init
    integer             ::  nspin
    integer             ::  ncell
    integer             ::  norb
    integer             ::  nsc

    real(8)             ::  Ebnd(2)
    integer             ::  estNe
    real(8)             ::  diag_acc
contains
    procedure       :: init => init_ham_init
end type

type parameters_TB_IO_H
    !parameters directly for the Hamiltonian
    !new parameters
    type(TB_hopping),allocatable    :: hop_io(:)   !hopping parameters
    type(TB_delta),allocatable      :: del_io(:)   !delta parameters
    type(TB_Jsd),allocatable        :: Jsd(:)   !delta parameters
    type(TBio_delta_onsite_scf),allocatable ::  del_scf_io(:)   !self-consistent delta on-site parameters
    type(TBio_defect),allocatable   ::  defect(:)    !defect onsite-terms in super-cell

    type(Htb_inp)       ::  hop
    type(Hdelta)        ::  del

    type(wann_dat)      ::  wann_io

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
    real(8)             ::  Ebnd(2)=[0.0d0,0.0d0]     !minimal and maximal energy values to consider in restricted eigensolver routines
    integer             ::  estNe=0                       !estimated number of eigenvalues in interval
    real(8)             ::  diag_acc=1d-12    ! accuracy of iterative eigenvalue solution (so far only fpm input)

    !selfconsistent delta parameters
    logical             ::  use_scf=.false.         !use selfconsistent delta
    logical             ::  scf_print=.false.       !print intermediate delta steps
    integer             ::  scf_loopmax=100         !maximal number of loop iterations converging delta
    real(8)             ::  scf_diffconv=1.0d-6     !convergence criterion for difference of delta sum
    real(8)             ::  scf_Ecut=-1.0d0         !energy cutoff for selfconsistent delta energy sum 
    integer             ::  scf_kgrid(3)=[10,10,1]  !kgrid for delta-selfconsistency cycle in reciprocal space

end type 

type parameters_TB_IO_EF
    real(8)     :: N_electrons=-1.0d0  !total number of electrons (negative value means use E_F) 
    real(8)     :: E_F_in=0.0d0        !fermi energy
    real(8)     :: kt=1.0d-3           !smearing factor in energy
end type

type parameters_TB_IO_DOS
    real(8)     :: E_ext(2)=[-1.0d0,1.0d0]          !minimum and maximum energy to plot in dos
    real(8)     :: dE=1.0d-2                        !energy binning size
    real(8)     :: sigma=1.0d-2                     !gauss smearing parameter for dos
    integer     :: kgrid(3)=[1,1,1]                 !number of k-points in each direction in case of k-space dos
    logical     :: print_kint=.false.               !print out the index of the currently considered k index 
    type(dos_bnd_io),allocatable    ::  bnd_io(:)   !io for local dos
    integer,allocatable :: bnd(:,:)                 !local dos bnd parameters (2,number local dos)
end type

type parameters_TB_IO_HIGHS
   integer      ::  N_highsym=-1               !Number of high symmetry points to be connected
   real(8),allocatable     ::  k_highs(:,:)    !(3,Nighsym): high symmetry kpoints that are to be connected
   integer,allocatable     ::  k_highs_frac(:) !(Nighsym): integer fractions to apply on input k_highs
   real(8)                 ::  aim_dist=1.0d-2 !aimed k-space distance between neighboring points along high symmetry line
end type

type parameters_TB_IO_OCC_MULT
   !parameters for multiple occupation calculation at different energies
   real(8)                 ::  dE=-1.0                      !energy step size to plot
   real(8)                 ::  E_ext(2)=[0.0d0,0.0d0]      !minimal and aimed maximal energy considered
   real(8)                 ::  kt=-1.0d0                    !smearing factor of fermi and derivative of fermi function
end type


type parameters_TB_Hsolve
    !might want to remove this and put it in the io_H or Hamiltonian-type
    !basic parameters
    integer         ::  nspin=1     !number of spins (1 or 2) for each orbital
    integer         ::  ncell=-1    !overall number of cells
    integer         ::  norb=-1     !number of orbitals in cell
    integer         ::  nsc=1       !2 if doubling for BdG superconductivity
    integer         ::  dimH=-1     !final size of Hamiltonian including all modifications
!    !externally set for compatibility with Energy routines
!    integer         ::  pos_ext(2)  !values to access outside parameters 
    !which method to use
    logical         ::  sparse=.False.
    integer         ::  i_diag=1  !different diagonalization methods
    logical         ::  rearrange=.false.  !rearrange Hamiltonian basis order to have same site c and c^+  next to each other
    real(8)         ::  diag_acc=1d-12    ! accuracy of iterative eigenvalue solution (so far only fpm input)

    !calculating only part of spectrum 
    real(8)         ::  Ebnd(2)=[0.0d0,0.0d0]     !minimal and maximal energy values to consider in restricted eigensolver routines
    integer         ::  estNe=0                    !estimated number of eigenvalues in interval  (0 correponds to dimH)

    contains
    procedure :: upd => upd_h_par

end type

type parameters_TB_IO_FLOW
    logical         ::  do_r=.False.
    logical         ::  dos_r=.False.
    logical         ::  spec_r=.False.
    logical         ::  fermi_r=.False.
    logical         ::  occ_r=.False.
    logical         ::  occ_mult_r=.False.

    logical         ::  read_solution_r=.True.
    logical         ::  write_solution_r=.False.

    logical         ::  do_k=.False.
    logical         ::  dos_k=.False.
    logical         ::  highs_k=.False.
    logical         ::  fermi_k=.False.
end type

type parameters_TB
    type(parameters_TB_IO_H)            ::  io_H
    type(parameters_TB_IO_EF)           ::  io_ef
    type(parameters_TB_IO_DOS)          ::  io_dos
    type(parameters_TB_IO_HIGHS)        ::  io_highs
    type(parameters_TB_IO_flow)         ::  flow
    type(parameters_TB_Hsolve)          ::  H
    type(parameters_TB_IO_OCC_MULT)     ::  io_occ_mult
    logical         ::  is_mag=.False. !Hamiltonian has spins
    logical         ::  is_sc=.False. !Hamiltonian is superconducting-> everything doubles to include creators and destructors
contains
    procedure :: init => init_parameters_TB
end type

contains

subroutine init_ham_init(this,io)
    class(parameters_ham_init),intent(out)  ::  this
    class(parameters_TB_IO_H),intent(in)    ::  io
    this%nspin   = io%nspin    
    this%ncell   = io%ncell   
    this%norb    = io%norb     
    this%nsc     = io%nsc     

    this%Ebnd    = io%Ebnd    
    this%estNe   = io%estNe    
    this%diag_acc= io%diag_acc
end subroutine

subroutine upd_h_par(this)
    class(parameters_TB_Hsolve) ::   this

    this%dimH=this%nsc*this%nspin*this%norb*this%ncell
end subroutine

subroutine init_parameters_TB(TB_params,lat)
    !subroutine which does operation to prepare the TB-hamiltonian 
    ! starting from the pure input using lattice properties etc.
    use m_derived_types, only: lattice
    class(parameters_TB),intent(inout)  :: TB_params
    type(lattice), intent(in)           :: lat

    integer :: i
    logical :: fexist

    !use superconducting version if any delta is set
    TB_params%is_sc=allocated(TB_params%io_H%del_io)
    !use magnetic version if any atom has both a magnetic moment and an orbital or if superconducting is set
    TB_params%is_mag=any(lat%cell%atomic%moment/=0.0d0.and.lat%cell%atomic%orbitals>0).or.TB_params%is_sc
    !set dimension of Hamiltonian parameters
    associate( par=>TB_params%io_H)
        if(TB_params%is_sc) par%nsc=2
        if(TB_params%is_mag) par%nspin=2
        par%ncell=lat%Ncell
        par%norb_at=lat%cell%atomic%orbitals
        par%norb_at_off=[(sum(par%norb_at(1:i-1)),i=1,size(par%norb_at))]
        par%norb=sum(par%norb_at)
        par%dimH=par%nsc*par%nspin*par%norb*par%ncell
        do i=1,size(par%hop_io)
            Call par%hop_io(i)%check(lat)
        enddo
        Call par%hop%set(par%hop_io,lat,par%nspin)
        inquire(file='delta_onsite.inp',exist=fexist)
        if(fexist)then
            Call par%del%read_file('delta_onsite.inp',lat)
        else
            if(allocated(par%del_io)) Call par%del%set(par%del_io,lat,par%norb_at_off)
        endif
        par%use_scf=allocated(par%del_scf_io)
        if(par%use_scf) Call par%del%set_scf(par%del_scf_io,lat,par%norb_at_off)
    end associate
    if(allocated(TB_params%io_dos%bnd_io))then
        do i=1,size(TB_params%io_dos%bnd_io)
            Call TB_params%io_dos%bnd_io(i)%check(lat)
        enddo
        Call dos_get_ind(TB_params%io_dos%bnd_io,lat,TB_params%io_H%nspin,TB_params%io_H%norb_at,TB_params%io_H%norb_at_off,TB_params%io_dos%bnd)
    endif

    !REMOVE-> move to H_io
    !set dimensions of the Hamiltonian
    if(TB_params%is_sc) TB_params%H%nsc=2
    if(TB_params%is_mag) TB_params%H%nspin=2
    TB_params%H%ncell=lat%Ncell
    TB_params%H%norb=sum(lat%cell%atomic%orbitals)
    Call TB_params%H%upd()

    TB_params%H%sparse=TB_params%io_H%sparse
    TB_params%H%i_diag=TB_params%io_H%i_diag
    TB_params%H%estNe=TB_params%io_H%estNe
    TB_params%H%Ebnd=TB_params%io_H%Ebnd
    TB_params%H%diag_acc=TB_params%io_H%diag_acc
end subroutine

end module
