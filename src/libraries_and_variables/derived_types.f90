module m_derived_types
use m_basic_types
use m_type_lattice

! parameters for printing in and out
type io_parameter
! Go in the spmstm program of Tobias
     logical :: io_spstmL=.false.
! plot the spin structure (or the order parameter structure)
     logical :: io_Xstruct=.false.
! plot the effective neighbouring field
     logical :: io_Beff=.false.
! plot the stochastic field
     logical :: io_Tfield=.false.
! frequency for writting the plotting data (magnetization density and so one)
     integer :: io_frequency
! plot the fourrier tranform of the spin structure (or the order parameter structure)
     logical :: io_fft_Xstruct=.false.
! plot the topological charge density distribution
     logical :: io_topo=.false.
! plot the emergent magnetic field
     logical :: io_topohall=.false.
! plot the warnings or not
     logical :: io_warning=.false.
! frequency of writting of the data in convergence.dat and EM.dat
     integer :: io_writing
! theta and phi distribution
     logical :: io_Angle_Distrib=.false.
! energy density distribution
     logical :: io_Energy_Distrib=.false.
! field density distribution
     logical :: io_Field_Distrib=.false.
! energy detail (keep more information about origin of energy terms) 
     logical :: io_Energy_detail=.false.
! write out energy contributions during the Efreq steps of the dynamics
     logical :: io_energy_cont=.false.
! force field density distribution
     logical :: io_Force=.false.
! Track singularities in a vector field
     logical :: io_tracker=.false.
end type io_parameter

! mpi variable
type parallel_variables
! for checkerboard parallelization
     logical :: io_ghost
     integer :: n_ghost
! for temperature parallelization
     logical :: i_separate,i_average,i_paratemp
end type parallel_variables

!!!!!!!!
! simulation parameters
!!!!!!!!

type simulation_parameters
    type(real_var) :: ktini=real_var(0.0d0,'initial temperature')
    type(real_var) :: ktfin=real_var(0.0d0,'final temperature')
    type(vec_var) :: H_ext=vec_var((/0.0d0,0.0d0,0.0d0/),'external B-field'),E_ext=vec_var((/0.0d0,0.0d0,0.0d0/),'external E-field')
end type simulation_parameters

end module m_derived_types
