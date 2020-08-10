module m_TB_types
implicit none
public

type parameters_TB_IO_H
    !parameters directly for the Hamiltonian
    real(kind=8), allocatable :: hopping(:,:,:)   ! 1: up or down, 2: orbital, 3: neighbor
    real(kind=8), allocatable :: onsite(:,:)     ! 1:up or down, 2: orbital
    logical :: is_magnetic=.false.
    logical :: is_sc=.false.
    integer :: nb_shell=-1
    integer :: nb_orbitals=-1
    integer :: nb_spin=1
    real(kind=8), allocatable :: Jsd(:)  !1: orbital
    complex(kind=8), allocatable :: delta(:)  !1: orbital   !super conductivity delta
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
    logical         ::  is_mag=.False. !Hamiltonian has spins
    logical         ::  is_sc=.False. !Hamiltonian is superconducting-> everything doubles to include creators and destructors
    !integer         :: dimH
    !real(8)         :: E_F
end type



end module
