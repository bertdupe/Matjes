module m_TB_types
implicit none
public
private upd_h_par
type parameters_TB_IO_H
    !parameters directly for the Hamiltonian
    real(kind=8), allocatable :: hopping(:,:,:)   ! 1: up or down, 2: orbital, 3: neighbor
    real(kind=8), allocatable :: onsite(:,:)     ! 1:up or down, 2: orbital
    integer :: nb_shell=-1
    integer :: nb_orbitals=-1
    integer :: nb_spin=1
    real(kind=8), allocatable :: Jsd(:)  !1: orbital
    complex(kind=8), allocatable :: delta(:)  !1: orbital   !super conductivity delta
    integer             ::  i_diag=1  !different diagonalization methods
    logical             ::  sparse=.false.  !do calculation sparse
end type 

type parameters_TB_IO_EF
    real(kind=8) :: N_electrons=-1.0d0  !total number of electrons (negative value means use E_F) 
    real(kind=8) :: E_F_in=0.0d0        !fermi energy
    real(kind=8) :: kt=1.0d-3           !smearing factor in energy
end type

type parameters_TB_IO_DOS
    real(8)                     :: E_ext(2)=[-1.0d0,1.0d0]      !minimum and maximum energy to plot in dos
    real(8)                     :: dE=1.0d-2                    !energy binning size
    real(8)                     :: sigma=1.0d-2                 !gauss smearing parameter for dos
end type

type parameters_TB_IO_HIGHS
   integer      ::  N_highsym=-1               !Number of high symmetry points to be connected
   real(8),allocatable     ::  k_highs(:,:)    !(3,Nighsym): high symmetry kpoints that are to be connected
   integer,allocatable     ::  k_highs_frac(:) !(Nighsym): integer fractions to apply on input k_highs
   real(8)                 ::  aim_dist=1.0d-2 !aimed k-space distance between neighboring points along high symmetry line
end type

type parameters_TB_Hsolve
    !basic parameters
    integer         ::  nspin=1     !number of spins (1 or 2) for each orbital
    integer         ::  ncell=-1    !overall number of cells
    integer         ::  norb=-1     !number of orbitals in cell
    integer         ::  nsc=1       !2 if doubling for BdG superconductivity
    !set through upd
    integer         ::  dimH=-1     !final size of Hamiltonian including all modifications
    integer         ::  nsite=1     !overall number of states in each cells
    !externally set for compatibility with Energy routines
    integer         ::  pos_ext(2)  !values to access outside parameters 
    !which method to use
    logical         ::  sparse=.False.
    integer         ::  i_diag=1  !different diagonalization methods

    contains
    procedure :: upd => upd_h_par

end type

type parameters_TB_IO_FLOW
    logical         ::  do_r=.False.
    logical         ::  dos_r=.False.
    logical         ::  spec_r=.False.
    logical         ::  fermi_r=.False.
    logical         ::  occ_r=.False.

    logical         ::  do_k=.False.
    logical         ::  dos_k=.False.
    logical         ::  highs_k=.False.
    logical         ::  fermi_k=.False.
end type

type parameters_TB
    type(parameters_TB_IO_H)      ::  io_H
    type(parameters_TB_IO_EF)     ::  io_ef
    type(parameters_TB_IO_DOS)    ::  io_dos
    type(parameters_TB_IO_HIGHS)  ::  io_highs
    type(parameters_TB_IO_flow)   ::  flow
    type(parameters_TB_Hsolve)     ::  H
    logical         ::  is_mag=.False. !Hamiltonian has spins
    logical         ::  is_sc=.False. !Hamiltonian is superconducting-> everything doubles to include creators and destructors
end type

contains
subroutine upd_h_par(this)
    class(parameters_TB_Hsolve) ::   this

    this%nsite=this%nsc*this%nspin*this%norb
    this%dimH=this%nsite*this%ncell
end subroutine

end module
