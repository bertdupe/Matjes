module m_derived_types
use m_basic_types

! unit cells
! the unit cell can be magnetic, ferroelectric or nothing and can also have transport
type cell
     real(kind=8), allocatable :: mom(:), pos(:,:)
     logical, allocatable :: i_mom(:)
end type cell

type lattice
     real(kind=8) :: areal(3,3),astar(3,3),alat(3)
     integer :: dim_lat(3),n_system
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

!Hamiltonian operator of the simulations
type total_Energy
     class(*), dimension(:,:,:), pointer :: Op
     integer :: nline,ncolumn
     integer :: neighbor
end type total_Energy

! Hamiltonian coefficients

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
! total Hamiltonian
     real(kind=8), allocatable, dimension(:,:,:) :: total
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
     logical :: i_ani=.true.
end type Coeff_Ham

!!!!!!!!
! simulation parameters
!!!!!!!!

type simulation_parameters
    type(real_var) :: kt=real_var(0.0d0,'temperature'),ktini=real_var(0.0d0,'initial temperature')
    type(real_var) :: ktfin=real_var(0.0d0,'final temperature')
    type(vec_var) :: H_ext=vec_var((/0.0d0,0.0d0,0.0d0/),'external B-field'),E_ext=vec_var((/0.0d0,0.0d0,0.0d0/),'external E-field')
end type simulation_parameters

!type site_ham_op
!     class(Ham_op), pointer, dimension(:) :: Op
!     integer :: nei
!end type site_ham_op





!Operators of the simulations
type operator_poly
    class (*), pointer, dimension(:,:) :: Op
    integer :: il,ic
    integer :: nline,ncolumn
end type operator_poly


!!!!
!Hamiltonian of the simulations
!!!!

       type H_loc
        real(kind=8), pointer :: H_loc(:,:)
       end type H_loc

!Hamiltonian of the simulations
       type real_local_Hamil
        real(kind=16), allocatable, dimension(:,:) :: coeff
       end type real_local_Hamil

       type im_local_Hamil
        complex(kind=16), allocatable, dimension(:,:) :: coeff
       end type im_local_Hamil


       end module m_derived_types
