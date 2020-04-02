!
!
!======================================================================
Program Matjes
use m_io_utils
use m_io_files_utils
use m_init_variables
use mtprng
use m_derived_types
use m_lattice
use m_write_spin
use m_createspinfile

! old code
      use m_vector, only : norm
      use m_average_MC
      use m_user_info
      use m_check_restart

#ifdef CPP_MPI
      use m_randperm
      use m_make_box
      use m_mpi_prop, only : irank,irank_working,isize,MPI_COMM,all_world,irank_box,MPI_COMM_MASTER,MPI_COMM_BOX
      use m_gather_reduce
#endif
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
type(cell) :: motif
! external parameter
type(simulation_parameters) :: ext_param
! internal variables
type(mtprng_state) :: state
! tag that defines the system
      integer :: n_system
      Integer :: N_cell
! the computation time
      real(kind=8) :: computation_time

  call welcome()

! initialize the random number generator
#ifdef CPP_MRG
#ifdef CPP_MPI
  call mtprng_init(irank+1,state)
#else
  call mtprng_init(1,state)
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

! prepare the different lattices for the simulations
!call setup_lattice(my_simu,all_lattices,io_simu)


! read the input and prepare the lattices, the Hamitlonian and all this mess
call setup_simu(my_simu,io_simu,all_lattices,motif,ext_param,state)

! number of cell in the simulation
N_cell=size(all_lattices%l_modes)
n_system=all_lattices%n_system

write(6,'(I6,a)') N_cell, ' unit cells'
write(6,'(a)') '-----------------------------------------------'
write(6,'(a)') ''
write(6,'(a)') '-----------------------------------------------'

call CreateSpinFile('Spinse_start.dat',all_lattices,motif)
call WriteSpinAndCorrFile('SpinSTM_start.dat',all_lattices)




!     Start main procedures:
!     *****************************************************************

!!!!!!!! part that does the tight-binging from a frozen spin configuration

!if (my_simu%name == 'tight-binding') then
!           call tightbinding(spin,shape_spin)
!endif
!!!!!!!! part of the parallel tempering

!if (my_simu%name == 'parallel-tempering') then
!            call parallel_tempering(i_biq,i_dip,i_DM,i_four,i_stone,ising,i_print_W,equi,overrel,sphere,underrel,cor_log,gra_log, &
!            &    spin,shape_spin,tableNN,shape_tableNN,masque,shape_masque,indexNN,shape_index,EA,n_system, &
!            &    n_sizerelax,T_auto,i_optTset,N_cell,print_relax,N_temp,T_relax_temp,kTfin,kTini,h_ext,cone,n_Tsteps, &
!            &    i_ghost,n_ghost,nRepProc,mag_lattice)

!            call cpu_time(computation_time)
!            write(*,*) 'computation time:',computation_time,'seconds'

!#ifdef CPP_MPI
!            if (irank_working.eq.0) Write(*,*) 'The program has reached the end.'
!            call MPI_FINALIZE(ierr)
!#else
!            Write(*,*) 'The program has reached the end.'
!#endif

!           stop
!endif

!---------------------------------
!  Part which does a normal MC with the metropolis algorithm
!---------------------------------

if (my_simu%name == 'metropolis') call MonteCarlo(all_lattices,motif,io_simu,ext_param)

!---------------------------------
!  Part which does the Spin dynamics
!    Loop for Spin dynamics
!---------------------------------

if (my_simu%name == 'magnet-dynamics') call spindynamics(all_lattices,motif,io_simu,ext_param)

!---------------------------------
!  Part which does Entropic Sampling
!---------------------------------
if (my_simu%name == 'entropics') call entropic(all_lattices,motif,io_simu,ext_param)


!---------------------------------
!  Part which does the GNEB
!---------------------------------
if (my_simu%name == 'GNEB') then
            write(6,'(a)') 'entering into the GNEB routine'
!!            call init_gneb()
             call GNEB(all_lattices,motif,io_simu,ext_param)
            !call set_gneb_defaults()
endif

!---------------------------------
!  Part which does the tight-binding simulations
!---------------------------------

if (my_simu%name == 'tight-binding') then
            write(6,'(a)') 'entering into the tight-binding routines'
             call tightbinding(all_lattices,motif,io_simu,ext_param)
endif


!if (my_simu%name == 'minimization') then
!            write(6,'(a)') 'entering into the minimization routine'
!
!            call CreateSpinFile('Spinse_start.dat',mag_lattice,mag_motif)
!
!            call minimize(i_biq,i_dm,i_four,i_dip,gra_log,gra_freq,EA, &
!              & spin,shape_spin,tableNN,shape_tableNN,masque,shape_masque,indexNN,shape_index,N_cell,h_ext,mag_lattice)
!endif ! montec, dynamic or i_gneb


!---------------------------------
!  Part which does the GNEB
!---------------------------------

!        if (my_simu%i_pimc) then
!            write(6,'(a)') 'entering into the path integral monte carlo routine'
!
!            call pimc(state,i_biq,i_dm,i_four,i_dip,EA,h_ext, &
!                    & spin,shape_spin,tableNN,shape_tableNN,masque,shape_masque,indexNN,shape_index,N_cell)

!        endif
!---------------------------------
!  Part which does the PIMC
!---------------------------------

call CreateSpinFile('Spinse_end.dat',all_lattices,motif)
call WriteSpinAndCorrFile('SpinSTM_end.dat',all_lattices)

call cpu_time(computation_time)
write(*,*) 'computation time:',computation_time,'seconds'
Write(*,*) 'The program has reached the end.'

    END
