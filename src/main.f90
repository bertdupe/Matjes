!
!
!======================================================================
Program Matjes
use m_io_utils
use m_io_files_utils
use m_init_variables
use m_derived_types
use m_lattice
use m_setup_simu
use m_write_spin
use m_createspinfile
use m_minimize
use m_get_random
use m_user_info
use m_H_public
use m_spindynamics
use m_montecarlo
use m_entropic
use m_GNEB, only: GNEB
use m_write_config, only: write_config
use m_rw_minimize, only: rw_minimize, min_input

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
! description of the unit cell
type(t_cell) :: motif
! external parameter
type(simulation_parameters) :: ext_param
! Hamiltonian used (extend to array with different basis + higher ranks)
class(t_H),allocatable      :: Ham_res(:), Ham_comb(:)
! tag that defines the system
      integer :: n_system
      Integer :: N_cell
! the computation time
      real(kind=8) :: computation_time
type(min_input)    :: io_min

  call welcome()

! initialize the random number generator
#ifdef CPP_MRG
#ifdef CPP_MPI
  call init_mtprng(irank+1)
#else
  call init_mtprng(1)
#endif
#else
  call init_rand_seed
#endif

io_param=open_file_read('input')

!! initialize the variables
! initialization of the type of simulation
call init_variables(my_simu)
! initialization of the lattice
call init_variables(all_lattices)
! initialization of the io_part
call init_variables(io_simu)

! get the simulation type
call get_parameter(io_param,'input',my_simu)
call close_file('input',io_param)

! read the input and prepare the lattices, the Hamitlonian and all this mess
call setup_simu(io_simu,all_lattices,motif,ext_param,Ham_res,Ham_comb)

! number of cell in the simulation
N_cell=size(all_lattices%ordpar%l_modes)
n_system=all_lattices%n_system

write(6,'(I6,a)') N_cell, ' unit cells'
write(6,'(a)') '-----------------------------------------------'
write(6,'(a)') ''
write(6,'(a)') '-----------------------------------------------'

!call CreateSpinFile('Spinse_start.dat',all_lattices,motif)
!call WriteSpinAndCorrFile('SpinSTM_start.dat',all_lattices)
Call write_config('start',all_lattices)




!     Start main procedures:
!     *****************************************************************

!!!!!!!! part of the parallel tempering

!if (my_simu%name == 'parallel-tempering') then
!            call parallel_tempering(i_biq,i_dip,i_DM,i_four,i_stone,ising,i_print_W,equi,overrel,sphere,underrel,cor_log,gra_log, &
!            &    spin,shape_spin,tableNN,shape_tableNN,masque,shape_masque,indexNN,shape_index,EA,n_system, &
!            &    n_sizerelax,T_auto,i_optTset,N_cell,print_relax,N_temp,T_relax_temp,kTfin,kTini,h_ext,cone,n_Tsteps, &
!            &    i_ghost,n_ghost,nRepProc,mag_lattice)

!            call cpu_time(computation_time)
!            write(*,*) 'computation time:',computation_time,'seconds'

!---------------------------------
!  Part which does a normal MC with the metropolis algorithm
!---------------------------------

if (my_simu%name == 'metropolis')then
    call MonteCarlo(all_lattices,io_simu,ext_param,Ham_comb)
endif

!---------------------------------
!  Part which does the Spin dynamics
!    Loop for Spin dynamics
!---------------------------------

if (my_simu%name == 'magnet-dynamics') call spindynamics(all_lattices,io_simu,ext_param,Ham_comb,Ham_res)

!---------------------------------
!  Part which does Entropic Sampling
!---------------------------------
if (my_simu%name == 'entropic')then
    STOP "put entropic back in"
    !call entropic(all_lattices,motif,io_simu,ext_param)
endif



!---------------------------------
!  Part which does the GNEB
!---------------------------------
if (my_simu%name == 'GNEB') then
            write(6,'(a)') 'entering into the GNEB routine'
            call GNEB(all_lattices,io_simu,ext_param,Ham_comb)
endif

!---------------------------------
!  Part which does the tight-binding simulations
!---------------------------------

if (my_simu%name == 'tight-binding') then
            write(6,'(a)') 'entering into the tight-binding routines'
             call tightbinding(all_lattices,io_simu)
endif

if (my_simu%name == 'minimization')then
    Call rw_minimize(io_min)
    call minimize(all_lattices,io_simu,io_min,Ham_comb)
endif

if (my_simu%name == 'minimize_infdamp')then
    Call rw_minimize(io_min)
    call minimize_infdamp(all_lattices,io_simu,io_min,Ham_comb)
endif
!---------------------------------
!  Part which does the molecular dynamics
!---------------------------------

if (my_simu%name == 'molecular_dynamics')then

     write(6,'(a)') 'entering into the molecular dynamics routines'
     call molecular_dynamics(all_lattices,io_simu)

endif

!---------------------------------
!---------------------------------

Call write_config('end',all_lattices)
!call CreateSpinFile('Spinse_end.dat',all_lattices,motif)
!call WriteSpinAndCorrFile('SpinSTM_end.dat',all_lattices)

call cpu_time(computation_time)
write(*,*) 'computation time:',computation_time,'seconds'
Write(*,*) 'The program has reached the end.'

    END