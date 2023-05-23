module m_derived_types
use m_basic_types
use m_type_lattice
use m_kgrid

! parameters for printing in and out
type io_parameter
! Go in the spmstm program of Tobias
     logical :: io_spstmL=.false.
! plot the spin structure (or the order parameter structure)
     logical :: io_Xstruct=.false.
! plot the effective neighbouring field
     logical :: io_Beff=.false.
! plot the forces
     logical :: io_Feff=.false.
! plot the stochastic field
     logical :: io_Tfield=.false.
! plot the fourrier tranform of the spin structure (or the order parameter structure)
     logical :: io_fft_Xstruct=.false.
! plot the Fourier transform of the Hamiltonian
     logical :: io_fft_Ham=.false.
! plot the topological charge density distribution
     logical :: io_topo=.false.
! plot the emergent magnetic field
     logical :: io_topohall=.false.
! plot the warnings or not
     logical :: io_warning=.false.
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
! frequency of writting of the data in convergence.dat and EM.dat
     integer :: io_writing=100
! frequency for writting the plotting data (magnetization density and so one)
     integer :: io_frequency=100
! calculate topological charged
    logical :: calc_topo=.true.
! calculate correlations
    logical :: calc_correlations=.true.
contains
    procedure :: bcast => io_parameter_bcast
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
    real(8) :: ktini=0.0d0
    real(8) :: ktfin=0.0d0
    real(8),dimension(3)    ::  H_ext=0.0d0,E_ext=0.0d0
    contains
    procedure :: bcast => simulation_parameters_bcast
    procedure :: set => simulation_parameters_set
end type simulation_parameters

contains

subroutine simulation_parameters_bcast(this,comm)
    use mpi_basic                
    class(simulation_parameters),intent(inout)      ::  this
    type(mpi_type),intent(in)                       ::  comm

#ifdef CPP_MPI
    integer     :: ierr
    Call MPI_Bcast(this%ktini, 1, MPI_REAL8, comm%mas, comm%com,ierr)
    Call MPI_Bcast(this%ktfin, 1, MPI_REAL8, comm%mas, comm%com,ierr)
    Call MPI_Bcast(this%H_ext, 3, MPI_REAL8, comm%mas, comm%com,ierr)
    Call MPI_Bcast(this%E_ext, 3, MPI_REAL8, comm%mas, comm%com,ierr)
#else
    continue
#endif
end subroutine

subroutine simulation_parameters_set(this,H,E,T)
    use m_constants, only : k_B
    class(simulation_parameters),intent(inout)      ::  this
    real(8),intent(in)  ::  H(3),E(3),T(2)

    this%H_ext=H
    this%E_ext=E
    this%ktini=T(1)*k_B
    this%ktfin=T(2)*k_B
end subroutine

subroutine io_parameter_bcast(this,comm)
    use mpi_basic                
    class(io_parameter),intent(inout)   ::  this
    type(mpi_type),intent(in)           ::  comm

#ifdef CPP_MPI
    integer     :: ierr
    Call MPI_Bcast(this%io_spstmL        , 1, MPI_LOGICAL, comm%mas, comm%com,ierr)
    Call MPI_Bcast(this%io_Xstruct       , 1, MPI_LOGICAL, comm%mas, comm%com,ierr)
    Call MPI_Bcast(this%io_Beff          , 1, MPI_LOGICAL, comm%mas, comm%com,ierr)
    Call MPI_Bcast(this%io_Tfield        , 1, MPI_LOGICAL, comm%mas, comm%com,ierr)
    Call MPI_Bcast(this%io_fft_Xstruct   , 1, MPI_LOGICAL, comm%mas, comm%com,ierr)
    Call MPI_Bcast(this%io_topo          , 1, MPI_LOGICAL, comm%mas, comm%com,ierr)
    Call MPI_Bcast(this%io_topohall      , 1, MPI_LOGICAL, comm%mas, comm%com,ierr)
    Call MPI_Bcast(this%io_warning       , 1, MPI_LOGICAL, comm%mas, comm%com,ierr)
    Call MPI_Bcast(this%io_Angle_Distrib , 1, MPI_LOGICAL, comm%mas, comm%com,ierr)
    Call MPI_Bcast(this%io_Energy_Distrib, 1, MPI_LOGICAL, comm%mas, comm%com,ierr)
    Call MPI_Bcast(this%io_Field_Distrib , 1, MPI_LOGICAL, comm%mas, comm%com,ierr)
    Call MPI_Bcast(this%io_Energy_detail , 1, MPI_LOGICAL, comm%mas, comm%com,ierr)
    Call MPI_Bcast(this%io_energy_cont   , 1, MPI_LOGICAL, comm%mas, comm%com,ierr)
    Call MPI_Bcast(this%io_Force         , 1, MPI_LOGICAL, comm%mas, comm%com,ierr)
    Call MPI_Bcast(this%io_tracker       , 1, MPI_LOGICAL, comm%mas, comm%com,ierr)
    Call MPI_Bcast(this%calc_topo        , 1, MPI_LOGICAL, comm%mas, comm%com,ierr)
    Call MPI_Bcast(this%io_writing       , 1, MPI_INTEGER, comm%mas, comm%com,ierr)
    Call MPI_Bcast(this%io_frequency     , 1, MPI_INTEGER, comm%mas, comm%com,ierr)
    !could be done more elegantly with custom MPI_type
#else
    continue
#endif
end subroutine

end module m_derived_types
