!
!
!======================================================================
      Program Matjes
      use m_write_spin
      use m_createspinfile
      use m_derived_types
      use m_Corre
      use m_constants
      use m_rw_lattice, only : motif,dim_lat,world,net,n_system
      use m_lattice, only : spin,tableNN,masque,indexNN
      use m_vector, only : norm
      use m_average_MC
      use m_topoplot
      use m_qorien
      use m_fft
      use mtprng
      use m_welcome
      use m_check_restart
      use m_parameters, only : n_Tsteps,n_sizerelax,i_qorien,CalTheta,Cor_log,Periodic_log,i_metropolis,dynamic,cone,gra_topo &
     & ,Total_MC_Steps,kt,n_thousand,T_auto,EA,i_gneb,gra_fft,CalEnergy,Gra_log,T_relax,kTfin,kTini,spstmL &
     & ,i_separate,i_average,i_ghost,i_paratemp,i_optTset,T_sweep,i_topohall,print_relax,i_minimization,gra_freq &
     & ,i_qorien,i_biq,i_dip,i_DM,i_four,i_stone,ising,i_print_W,equi,overrel,sphere,underrel,h_ext,N_temp,T_relax_temp,n_ghost &
     & ,nRepProc,i_entropic
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
! types of the simulation
      type(type_simu) :: my_simu
! for the random number generator
      type(mtprng_state) :: state
      integer :: alloc_check
! slope index
      Integer :: i_MC, n_MC, n_kt,i_x,i_y,i_z,i_m,k
! size of the different table
      integer :: shape_index(2),shape_spin(5),shape_tableNN(6),shape_masque(4)
!     slope index for corre, later to write files
      Integer :: i_lat,LenStr
      Integer :: j_MC,i,N_cell
!     name of the file
      character(len=50) :: fname,fname2,fname3
      character(len=26) :: toto
      logical :: exists,godsky,i_restart
      real(kind=8) :: l(3)
      Integer :: iii
      real(kind=8) :: kti
! sweep temperature configuration or have a geometric Temperature set
      logical :: i_sweep,i_geometric
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

#ifdef CPP_MPI
      if (irank.eq.0) call welcome()
#else
      call welcome()
#endif
#ifdef CPP_MRG
#ifdef CPP_MPI
      call mtprng_init(irank+1,state)
#else
      call mtprng_init(1,state)
#endif
#else
      call init_rand_seed
#endif
! read the input
      call setup_simu(my_simu,io_simu,state)

! number of cell in the simulation
      N_cell=count(masque(1,:,:,:).ne.0)
#ifdef CPP_MPI
      if (irank_working.eq.0) write(6,'(I6,a)') N_cell, ' unit cells'
#else
      write(6,'(I6,a)') N_cell, ' unit cells'
#endif

! ******************************
! Prepare the addition or removal of skyrmions.

    inquire (file='rwskyrmion.in',exist=godsky)
    if (godsky) then

      if (i_print_W) then
           WRITE(6,'(/,a)') 'WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!'
           WRITE(6,'(a)') 'You are creating artificially topologically protected spin structure'
           WRITE(6,'(a)') 'Is that really what you want? If not, please have'
           WRITE(6,'(a,/)') 'a look at lines 373 and 374'
      endif

        call rw_skyrmion()

    end if



! *****************************
! *****************************
! checking for spin dynamics

!    dynamic=.False.
!    inquire (file='dyna.in',exist=dynamic)
!    if (dynamic) then
#ifdef CPP_MPI
!        if (irank_working.eq.0) write(6,'(a)') 'Spin dynamics selected'
#else
!        write(6,'(a)') 'Spin dynamics selected'
#endif
!        i_metropolis=.False.
!    endif
! check for the use of teh GNEB
!    inquire (file='GNEB.in',exist=i_gneb)
!    if (i_gneb) then
!     write(6,'(a)') 'you are using GNEB method'
!     write(6,'(a)') 'you will not use MC or spin dynamics'
!     i_metropolis=.False.
!     dynamic=.False.
!    endif
! check for the minimization routine
!    inquire (file='minimization.in',exist=i_minimization)
!    if (i_minimization) then
!     close(7)
!     write(6,'(a)') 'you are using only the minimization routine of one spin structure'
!     write(6,'(a)') 'you will not use MC or spin dynamics or GNEB'
!     write(6,'(a)') 'This algorithms works only serial'
!     i_metropolis=.False.
!     dynamic=.False.
!     i_gneb=.False.
!    endif

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


!!!!!!!!!!
!!! initialize the variable for the size of the different table
    shape_index=shape(indexNN)
    shape_spin=shape(spin)
    shape_tableNN=shape(tableNN)
    shape_masque=shape(masque)

! ***********************
! Here you can choose to add or remove a skyrmion.

        if (godsky) call kornevien(kt)

!!!!!!!! part of the parallel tempering

        if (my_simu%i_paratemp) then
            call parallel_tempering(i_biq,i_dip,i_DM,i_four,i_stone,ising,i_print_W,equi,overrel,sphere,underrel,cor_log,gra_log, &
            &    spin,shape_spin,tableNN,shape_tableNN,masque,shape_masque,indexNN,shape_index,EA,n_system, &
            &    n_sizerelax,T_auto,i_optTset,N_cell,print_relax,N_temp,T_relax_temp,kTfin,kTini,h_ext,cone,state,n_Tsteps, &
            &    i_ghost,n_ghost,nRepProc,size(world))

            call cpu_time(computation_time)
            write(*,*) 'computation time:',computation_time,'seconds'

#ifdef CPP_MPI
            if (irank_working.eq.0) Write(*,*) 'The program has reached the end.'
            call MPI_FINALIZE(ierr)
#else
            Write(*,*) 'The program has reached the end.'
#endif

            stop
        endif

!---------------------------------
!  Part which does a normal MC with the metropolis algorithm
!---------------------------------
        if (my_simu%i_metropolis) call MonteCarlo(N_cell,n_system,state,size(world),n_thousand,n_sizerelax, &
            &    spin,shape_spin,tableNN,shape_tableNN,masque,shape_masque,indexNN,shape_index, &
            &    i_qorien,i_biq,i_dip,i_DM,i_four,i_stone,ising,i_print_W,equi,overrel,sphere,underrel,i_ghost, &
            &    Periodic_log,gra_topo,CalEnergy,CalTheta,Gra_log,spstmL,gra_fft, &
            &    i_separate,i_average,i_topohall,print_relax,Cor_log, &
            &    n_Tsteps,cone,Total_MC_Steps,T_auto,EA,T_relax,kTfin,kTini,h_ext, &
            &    n_ghost,nRepProc)
        
!---------------------------------
!  Part which does the Spin dynamics
!    Loop for Spin dynamics
!---------------------------------
        if (my_simu%i_dynamic) call spindynamics(state,spin,shape_spin,tableNN,shape_tableNN,masque,shape_masque,indexNN,shape_index,world, &
      &  N_cell, &
      &  i_biq,i_dm,i_four,i_dip,gra_topo,gra_log,gra_freq, &
      &  ktini,ktfin,EA,h_ext)

!---------------------------------
!  Part which does Entropic Sampling
!---------------------------------
        if (my_simu%i_entropic) call entropic(N_cell,n_system,state,size(world), &
            &    spin,shape_spin,tableNN,shape_tableNN,masque,shape_masque,indexNN,shape_index, &
            &    i_biq,i_dip,i_DM,i_four,i_stone,ising,i_print_W,equi,sphere,overrel,underrel,i_ghost, &
            &    Periodic_log,Gra_log,spstmL,gra_fft,cone,EA,T_relax,h_ext,n_ghost, &
            &    kTfin,kTini,n_Tsteps)


!---------------------------------
!  Part which does the GNEB
!---------------------------------
        if (my_simu%i_gneb) then
            write(6,'(a)') 'entering into the GNEB routine'
!            call init_gneb()
            call GNEB(state,i_biq,i_dm,i_four,i_dip,EA,h_ext, &
                    & indexNN,shape_index,shape_spin,tableNN,shape_tableNN,masque,shape_masque,N_cell)
            !call set_gneb_defaults()
        endif

        if (my_simu%i_minimization) then
            write(6,'(a)') 'entering into the minimization routine'

            call CreateSpinFile('Spinse_start.dat',spin,shape_spin)

            call minimize(i_biq,i_dm,i_four,i_dip,gra_log,gra_freq,EA, &
              & spin,shape_spin,tableNN,shape_tableNN,masque,shape_masque,indexNN,shape_index,N_cell,h_ext)
        endif ! montec, dynamic or i_gneb


!---------------------------------
!  Part which does the GNEB
!---------------------------------

        if (my_simu%i_pimc) then
            write(6,'(a)') 'entering into the path integral monte carlo routine'
!            call init_gneb()
            call pimc(state,i_biq,i_dm,i_four,i_dip,EA,h_ext, &
                    & spin,shape_spin,tableNN,shape_tableNN,masque,shape_masque,indexNN,shape_index,N_cell)
            !call set_gneb_defaults()
        endif
!---------------------------------
!  Part which does the PIMC
!---------------------------------

#ifndef CPP_MPI
! write the last spin structure
        call WriteSpinAndCorrFile('SpinSTM_end.dat',spin,shape_spin)

        call CreateSpinFile('Spinse_end.dat',spin,shape_spin)
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
