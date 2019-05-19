!
!
!======================================================================
Program Matjes
use m_io_utils
use m_io_files_utils
use m_init_variables
use mtprng
use m_derived_types

! old code
      use m_write_spin
      use m_createspinfile
      use m_Corre
      use m_constants
      use m_vector, only : norm
      use m_average_MC
      use m_topoplot
      use m_qorien
      use m_fft
      use m_user_info
      use m_check_restart
      use m_parameters, only : n_Tsteps,n_sizerelax,i_qorien,CalTheta,Cor_log,cone,gra_topo &
     & ,Total_MC_Steps,n_thousand,T_auto,gra_fft,CalEnergy,Gra_log,T_relax,spstmL &
     & ,i_separate,i_average,i_ghost,i_optTset,i_topohall,print_relax,gra_freq &
     & ,i_qorien,i_print_W,equi,overrel,sphere,underrel,N_temp,T_relax_temp,n_ghost &
     & ,nRepProc
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
Integer :: nlattice
type(lattice) :: all_lattices
! unit parameters for the file
integer :: io_param
! description of the unit cell
type(cell) :: mag_motif
! external parameter
type(simulation_parameters) :: ext_param
! internal variables
type(mtprng_state) :: state
! definition of the motif and the lattice
type(lattice) :: mag_lattice
! tag that defines the system
      integer :: n_system
! size of the different table
      integer :: shape_index(2),shape_spin(5),shape_tableNN(6),shape_masque(4)
      Integer :: N_cell
! the computation time
      real(kind=8) :: computation_time

#ifdef CPP_MPI
      integer :: ierr
      logical :: i_check

      include 'mpif.h'

      CALL MPI_INIT(ierr)
      call mpi_comm_group(MPI_COMM_WORLD,all_world,ierr)
      call mpi_comm_create(MPI_COMM_WORLD,all_world,MPI_COMM,ierr)
      call MPI_COMM_RANK(MPI_COMM, irank, ierr)
      CALL MPI_COMM_SIZE(MPI_COMM,isize,ierr)
      if (isize.lt.12) write(6,'(a,2x,I6,2x,a)')'MPI task',irank,'has started...'
      if (irank.eq.0) call welcome(isize)
#endif

! welcome message
#ifdef CPP_MPI
  if (irank.eq.0) call welcome()
#else
  call welcome()
#endif

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

call get_parameter(io_param,'input','nlattice',nlattice)

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


! read the input
call setup_simu(my_simu,io_simu,all_lattices,mag_motif,ext_param,state)

! number of cell in the simulation
N_cell=size(mag_lattice%l_modes)
n_system=mag_lattice%n_system

#ifdef CPP_MPI
      if (irank_working.eq.0) write(6,'(I6,a)') N_cell, ' unit cells'
#else
      write(6,'(I6,a)') N_cell, ' unit cells'
#endif

#ifdef CPP_MPI
    if (irank_working.eq.0) then
#endif
     write(6,'(a)') '-----------------------------------------------'
     write(6,'(a)') ''
     write(6,'(a)') '-----------------------------------------------'
#ifdef CPP_MPI
    endif
#endif






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

if (my_simu%name == 'metropolis') call MonteCarlo(all_lattices,mag_motif,io_simu,gra_topo,ext_param)
!   call MonteCarlo(N_cell,n_thousand,n_sizerelax, &
!            &    spin,shape_spin,tableNN,shape_tableNN,masque,shape_masque,indexNN,shape_index, &
!            &    i_qorien,i_biq,i_dip,i_DM,i_four,i_stone,ising,i_print_W,equi,overrel,sphere,underrel,i_ghost, &
!            &    gra_topo,CalEnergy,CalTheta,Gra_log,spstmL,gra_fft, &
!            &    i_separate,i_average,i_topohall,print_relax,Cor_log, &
!            &    n_Tsteps,cone,Total_MC_Steps,T_auto,EA,T_relax,kTfin,kTini,h_ext, &
!            &    n_ghost,nRepProc,mag_lattice)

!---------------------------------
!  Part which does the Spin dynamics
!    Loop for Spin dynamics
!---------------------------------

if (my_simu%name == 'magnet-dynamics') call spindynamics(all_lattices,mag_motif,io_simu,gra_topo,ext_param)

!---------------------------------
!  Part which does Entropic Sampling
!---------------------------------
!if (my_simu%name == 'entropics') call entropic(N_cell,n_system, &
!            &    spin,shape_spin,tableNN,shape_tableNN,masque,shape_masque,indexNN,shape_index, &
!            &    i_biq,i_dip,i_DM,i_four,i_stone,ising,i_print_W,equi,sphere,overrel,underrel,i_ghost, &
!            &    spstmL,gra_fft,cone,EA,T_relax,h_ext,n_ghost, &
!            &    kTfin,kTini,n_Tsteps,mag_lattice)


!---------------------------------
!  Part which does the GNEB
!---------------------------------
!if (my_simu%name == 'GNEB') then
!            write(6,'(a)') 'entering into the GNEB routine'
!!            call init_gneb()
!            call GNEB(i_biq,i_dm,i_four,i_dip,EA,h_ext,mag_lattice, &
!                    & indexNN,shape_index,shape_spin,tableNN,shape_tableNN,masque,shape_masque,N_cell)
!            !call set_gneb_defaults()
!endif

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

#ifndef CPP_MPI
! write the last spin structure
        call WriteSpinAndCorrFile(all_lattices)

        call CreateSpinFile(all_lattices,mag_motif)
#endif

        call cpu_time(computation_time)
        write(*,*) 'computation time:',computation_time,'seconds'

#ifdef CPP_MPI
        if (irank_working.eq.0) Write(*,*) 'The program has reached the end.'
        call MPI_FINALIZE(ierr)
#else
        Write(*,*) 'The program has reached the end.'
#endif
    END
