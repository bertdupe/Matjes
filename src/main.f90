!
!
!======================================================================
Program Matjes
use m_io_utils
use m_io_files_utils
use m_derived_types
use m_setup_simu
use m_write_spin
use m_createspinfile
use m_minimize
use m_llg_diag
use m_user_info
use m_H_public
use m_parallel_tempering
use m_spindynamics
use m_molecular_dynamics
use m_montecarlo
use m_entropic
use m_GNEB, only: GNEB
use m_write_config, only: write_config
use m_random_init, only: random_init
use m_hamiltonian_collection, only: hamiltonian

use m_mpi_start_end  !also includes mpi_basic
use m_fftw3,only: fftw_init

Implicit None

! variable for the simulation
!----------------------------------------
! io_variables
    type(io_parameter) :: io_simu
! types of the simulation (it defines what type of simulation will be run)
    type(bool_var) :: my_simu
! number of lattices that should be used + variables that contains all lattices
    type(lattice) :: all_lattices
! unit parameters for the file
    integer :: io_param
! external parameter
    type(simulation_parameters) :: ext_param
! Hamiltonian used (extend to array with different basis + higher ranks)
    class(t_H),allocatable      :: Ham_res(:), Ham_comb(:)
    type(hamiltonian)           :: H_res,H_comb !newer Hamiltonian class including the FFT parts of the Hamiltonian as well
! the computation time
    real(kind=8) :: computation_time

    Call init_MPI(mpi_world)
    Call fftw_init()

    call welcome()

! initialize the random number generator
    Call random_init()

    if(mpi_world%ismas)then
        ! get the simulation type
        io_param=open_file_read('input')
        call get_parameter(io_param,'input',my_simu)
        call close_file('input',io_param)
        
        ! read the input and prepare the lattices, the Hamitlonian and all this mess
        call setup_simu(io_simu,all_lattices,ext_param,Ham_res,Ham_comb,H_res,H_comb)
        
        Call write_config('start',all_lattices) !write the initial configuration of all considered states
    endif
    
    !     Start main procedures:
    !     *****************************************************************
    Call my_simu%bcast(mpi_world)
 

    !---------------------------------
    !  Part which does the parallel tempering
    !---------------------------------
    if (my_simu%name == 'parallel-tempering')then
        Call parallel_tempering(all_lattices,io_simu,ext_param,H_comb,mpi_world)
    endif

    
    !---------------------------------
    !  Part which does a normal MC with the metropolis algorithm
    !---------------------------------
    if (my_simu%name == 'metropolis')then
        call MonteCarlo(all_lattices,io_simu,ext_param,H_res,mpi_world)
    endif

    !---------------------------------
    !  Part which does spin-dynamics
    !---------------------------------
    if (my_simu%name == 'magnet-dynamics')then
        call spindynamics(all_lattices,io_simu,ext_param,H_comb,H_res,mpi_world)
    endif


    !---------------------------------
    !  Part which does the GNEB
    !---------------------------------
    if (my_simu%name == 'GNEB') then
        call GNEB(all_lattices,io_simu,H_comb,mpi_world)  !no parallelization implemented
    endif

    !---------------------------------
    !  minimizations
    !---------------------------------
    if (my_simu%name == 'minimization')then
        call minimize(all_lattices,io_simu,H_comb,mpi_world)  !no parallelization implemented
    endif
    
    if (my_simu%name == 'minimize_infdamp')then
        call minimize_infdamp(all_lattices,io_simu,H_comb,mpi_world)  !no parallelization implemented
    endif



    if(mpi_world%ismas)then
        !---------------------------------
        !  Part which does the Spin dynamics
        !    Loop for Spin dynamics
        !---------------------------------
        !---------------------------------
        !  Part which does Entropic Sampling
        !---------------------------------
        if (my_simu%name == 'entropic')then
            STOP "entropic currently not implemented"
            !call entropic(all_lattices,all_lattices%cell,io_simu)
        endif
        
        !---------------------------------
        !  Part which does the tight-binding simulations
        !---------------------------------
        
        if (my_simu%name == 'tight-binding') then
            write(6,'(a)') 'entering into the tight-binding routines'
            call tightbinding(all_lattices,io_simu)
        endif
        
        !---------------------------------
        !  diag llg
        !---------------------------------
        if (my_simu%name == 'llg_diag')then
            call diag_llg(all_lattices,io_simu,H_comb)
        endif
        
        !---------------------------------
        !  Part which does the molecular dynamics
        !---------------------------------
        
        if (my_simu%name == 'molecular_dynamics')then
            write(6,'(a)') 'entering into the molecular dynamics routines'
            call molecular_dynamics(all_lattices,io_simu,ext_param,Ham_comb)
        endif
        
        !---------------------------------
        !---------------------------------
        
        Call write_config('end',all_lattices)
    
        call cpu_time(computation_time)
        write(*,*) 'computation time:',computation_time,'seconds'
        Write(*,*) 'The program has reached the end.'
    endif
    Call end_MPI()
END Program
