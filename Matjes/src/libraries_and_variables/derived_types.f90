module m_derived_types
use m_basic_types

! unit cells
! the unit cell can be magnetic, ferroelectric or nothing and can also have transport
type cell
     type(atom), allocatable :: atomic(:)
     real(kind=8), allocatable :: mom(:), pos(:,:)
     logical, allocatable :: i_mom(:)
end type cell

! variable that defines the lattice
type lattice
     real(kind=8) :: areal(3,3),astar(3,3),alat(3)
     integer :: dim_lat(3),n_system,dim_mode
     integer, allocatable :: world(:)
     logical :: boundary(3)
! Table of pointer
     type(vec_point), allocatable :: l_modes(:,:,:,:)
end type lattice

! parameters for printing in and out
type io_parameter
     logical :: io_dispersion,io_qorien
! Go in the spmstm program of Tobias
     logical :: io_spstmL, io_spstmonly
! plot the spin structure (or the order parameter structure)
     logical :: io_Xstruct
! plot the stochastic field
     logical :: io_Tfield
! frequency for writting the plotting data (magnetization density and so one)
     integer :: io_frequency
! plot the fourrier tranform of the spin structure (or the order parameter structure)
     logical :: io_fft_Xstruct
! plot the topological charge density distribution
     logical :: io_topo
! plot the emergent magnetic field
     logical :: io_topohall
! plot the warnings or not
     logical :: io_warning
! frequency of writting of the data in convergence.dat and EM.dat
     integer :: io_writing
! theta and phi distribution
     logical :: io_Angle_Distrib
! energy density distribution
     logical :: io_Energy_Distrib
! field density distribution
     logical :: io_Field_Distrib
end type io_parameter

! mpi variable
type parallel_variables
! for checkerboard parallelization
     logical :: io_ghost
     integer :: n_ghost
! for temperature parallelization
     logical :: i_separate,i_average,i_paratemp
end type parallel_variables

! Hamiltonian variables
type site_Ham
     real(kind=8), dimension(:,:), allocatable :: H
end type site_Ham

type shell_Ham
     type(site_Ham), dimension(:), allocatable :: atom
end type

type point_shell_Operator
     type(Op_real), allocatable, dimension(:) :: shell
end type point_shell_Operator

type point_shell_mode
     type(vec_point), allocatable, dimension(:) :: shell
end type point_shell_mode

!!!!!!!!
! Hamiltonian type
!!!!!!!!

type coeff_ham_inter_spec
     real(kind=8) :: c_ham=-1.0d0
     character(len=30) :: name=''
     logical :: i_exist=.false.
     type(site_Ham), allocatable, dimension(:) :: ham
end type coeff_ham_inter_spec

type coeff_ham_inter_spec_pointer
     character(len=30) :: name=''
     type(Op_real), allocatable :: ham(:)
end type coeff_ham_inter_spec_pointer

! Hamiltonian coefficients
! to be deleted
type Coeff_Ham
     real(kind=8) :: c_Ji=-1.0d0
     real(kind=8) :: c_DM=-1.0d0
     real(kind=8) :: c_JB=-1.0d0
     real(kind=8) :: c_Ki=-1.0d0
     real(kind=8) :: c_ani=1.0d0
! Exchange interaction
     real(kind=8), allocatable, dimension(:,:,:) :: exchange
! DMI interaction
     real(kind=8), allocatable, dimension(:,:,:) :: DMI
! magnetocrystalline interaction
     real(kind=8), allocatable, dimension(:,:,:) :: ani
! magnetocrystalline interaction
     real(kind=8), allocatable, dimension(:,:) :: Zeeman
! total Hamiltonian
     type(shell_Ham), allocatable, dimension(:) :: total_shell
! stoner parameter
     real(kind=8) :: Ist=0.0d0
! biquadratic interaction
     real(kind=8) :: Biq=0.0d0
! 4-spin interaction
     real(kind=8) :: fours=0.0d0
! presence or absence of interactions
     logical :: i_DM=.false.
     logical :: i_four=.false.
     logical :: i_biq=.false.
     logical :: i_dip=.false.
     logical :: i_exch=.false.
     logical :: i_ani=.false.
end type Coeff_Ham

!!!!!!!!
! simulation parameters
!!!!!!!!

type simulation_parameters
    type(real_var) :: ktini=real_var(0.0d0,'initial temperature')
    type(real_var) :: ktfin=real_var(0.0d0,'final temperature')
    type(vec_var) :: H_ext=vec_var((/0.0d0,0.0d0,0.0d0/),'external B-field'),E_ext=vec_var((/0.0d0,0.0d0,0.0d0/),'external E-field')
end type simulation_parameters

!!!!!!!!!
! operator type
!!!!!!!!!
type operator_real
    type(Op_real), allocatable, dimension(:,:) :: value
    integer :: nline,ncolumn
    integer, allocatable :: line(:,:)
end type operator_real

!!!!!!!!!
! order parameter type
!!!!!!!!!
type order_parameter
    character(len=30) :: name
    integer :: start,end
end type

end module m_derived_types
