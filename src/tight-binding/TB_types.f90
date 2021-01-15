module m_TB_types
use m_input_H_types
use m_types_tb_h_inp
use m_ham_arrange
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
    type(TB_delta),allocatable      :: del(:)   !delta parameters
    type(TB_Jsd),allocatable        :: Jsd(:)   !delta parameters

    type(Htb_inp)       ::  hop

    integer             ::  nspin=1         !number of spins (1 or 2) for each orbital
    integer             ::  ncell=-1        !overall number of cells
    integer             ::  norb=-1         !number of orbitals in cell
    integer             ::  nsc=1           !2 if doubling for BdG superconductivity
    integer             ::  dimH=-1         !final size of Hamiltonian including all modifications
    integer,allocatable ::  norb_at(:)      !number of orbitals at each atom
    integer,allocatable ::  norb_at_off(:)  !offset of orbitals at each atom

    integer             ::  i_diag=1  !different diagonalization methods
    logical             ::  sparse=.false.  !do calculation sparse
    logical             ::  rearrange=.false.  !rearrange Hamiltonian basis order to have same site c and c^+  next to each other

    real(8)             ::  Ebnd(2)=[0.0d0,0.0d0]     !minimal and maximal energy values to consider in restricted eigensolver routines
    integer             ::  estNe=0                       !estimated number of eigenvalues in interval
    real(8)             ::  diag_acc=1d-12    ! accuracy of iterative eigenvalue solution (so far only fpm input)
end type 

type parameters_TB_IO_EF
    real(8)     :: N_electrons=-1.0d0  !total number of electrons (negative value means use E_F) 
    real(8)     :: E_F_in=0.0d0        !fermi energy
    real(8)     :: kt=1.0d-3           !smearing factor in energy
end type

type parameters_TB_IO_DOS
    real(8)     :: E_ext(2)=[-1.0d0,1.0d0]      !minimum and maximum energy to plot in dos
    real(8)     :: dE=1.0d-2                    !energy binning size
    real(8)     :: sigma=1.0d-2                 !gauss smearing parameter for dos
    integer     :: kgrid(3)=[1,1,1]    !number of k-points in each direction in case of k-space dos
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



end module
