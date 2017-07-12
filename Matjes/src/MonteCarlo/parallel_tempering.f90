      subroutine parallel_tempering(i_biq,i_dip,i_DM,i_four,i_stone,ising,i_print_W,equi,overrel,sphere,underrel,cor_log,gra_log, &
            &    spin,shape_spin,tableNN,shape_tableNN,masque,shape_masque,indexNN,shape_index,EA,n_system, &
            &    n_sizerelax,T_auto,i_optTset,N_cell,print_relax,N_temp,T_relax_temp,kTfin,kTini,h_ext,coni,state,n_Tsteps, &
            &    i_ghost,n_ghost,nRepProc,world)
      use mtprng
      use m_topocharge_all
      use m_set_temp
      use m_average_MC
      use m_constants, only : k_b
      use m_vector, only : norm
      use m_paratemp
      use m_store_relaxation
      use m_check_restart
      use m_createspinfile
#ifdef CPP_MPI
      use m_mpi_prop, only : irank_working,isize,MPI_COMM,all_world,irank_box,MPI_COMM_MASTER,MPI_COMM_BOX
      use m_mpi
#endif
      implicit none
      interface
        subroutine MCsteps(state,E_total,E_decompose,Magnetization,kt,acc,rate,tries,cone,n_system, &
            & spin,shape_spin,tableNN,shape_tableNN,masque,shape_masque,indexNN,shape_index,h_ext,EA, &
            & i_biq,i_dip,i_DM,i_four,i_stone,ising,i_print_W,equi,overrel,sphere,underrel,n_world)
          use mtprng
          integer, intent(in) :: shape_index(2),shape_spin(5),shape_tableNN(6),shape_masque(4),n_system
          integer, intent(in) :: tableNN(shape_tableNN(1),shape_tableNN(2),shape_tableNN(3),shape_tableNN(4),shape_tableNN(5),shape_tableNN(6))
          integer, intent(in) :: masque(shape_masque(1),shape_masque(2),shape_masque(3),shape_masque(4))
          integer, intent(in) :: indexNN(shape_index(1),shape_index(2)),n_world
          logical, intent(in) :: i_biq,i_dip,i_DM,i_four,i_stone,ising,i_print_W,equi,overrel,sphere,underrel
          real(kind=8), intent(in) :: kt,h_ext(3),EA(3)
          real(kind=8), intent(inout) :: spin(shape_spin(1),shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5))
          type(mtprng_state), intent(inout) :: state
          real(kind=8), intent(inout) :: E_total,Magnetization(3),E_decompose(8),acc,rate,cone,tries
        end subroutine
      end interface
      logical, intent(in) :: i_biq,i_dip,i_DM,i_four,i_stone,ising,i_print_W,equi,overrel,sphere,underrel,cor_log,i_ghost,gra_log
      integer, intent(in) :: N_cell,n_sizerelax,N_temp,T_relax_temp,n_Tsteps,n_system,T_auto,n_ghost,nRepProc,world
      integer, intent(in) :: shape_spin(5),shape_tableNN(6),shape_masque(4),shape_index(2)
      integer, intent(in) :: tableNN(shape_tableNN(1),shape_tableNN(2),shape_tableNN(3),shape_tableNN(4),shape_tableNN(5),shape_tableNN(6))
      integer, intent(in) :: masque(shape_masque(1),shape_masque(2),shape_masque(3),shape_masque(4))
      integer, intent(in) :: indexNN(shape_index(1),shape_index(2))
      real(kind=8), intent(in) :: kTfin,kTini,h_ext(3),coni,EA(3)
      real(kind=8), intent(inout) :: spin(shape_spin(1),shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5))
      type(mtprng_state), intent(inout) :: state
      logical, intent(in) :: i_optTset,print_relax
! internal variable
! parallel tempering label. can take 3 values
! 0: no label; +1: up; -1: down
      real(kind=8), allocatable :: label_world(:),nup(:),ndown(:),success_PT(:)
! convergence parameters
      real(kind=8) :: rate,cone,Nsuccess,tries
! all the temperatures
      real(kind=8),allocatable :: kT_all(:),kT_updated(:)  ! contains the temperatures T/k_b: kt_all contains the temperatures from the smalles to the largest
      integer,allocatable :: image_temp(:)  ! contains the POSITION of temperatures
! slope of temperature sets
      integer :: j_optset
! slope of the MC
      integer :: i_MC,i_thousand,i_relax,i_pos
      real(kind=8) :: pos
! in case if restart
      integer :: restart_MC_steps,restart_n_thousand
! size of all the tables
      integer :: size_table, size_M_replica
! internal variable
! variables that being followed during the simulation
      real(kind=8) :: qeulerp,qeulerm,vortex(3),magnetization(3),E_total
! keep the information of the energies of each replicas before each swap
      real(kind=8), allocatable :: E_temp(:)
! contribution of the different energy parts
      real(kind=8) :: E_decompose(8)
! slope for the matrices
      integer :: i_x,i_y,i_z,i_m
! considered temperature
      real(kind=8) :: kT
! autocorrelation for the parallel tempering and number of relaxation steps
      integer :: autocor_steps,relaxation_steps,total_relax_steps,freq_rw,i_rw
! slope for the number of images and the temperature
      integer :: i_image,i_temp
! restart of the parallel tempering
      logical :: i_reread_paratemp,i_rewrite_paratemp
! restart of the measured data
      logical :: i_reread_data,i_rewrite_data
! restart of of the replicas
      logical :: i_reread_replicas,i_rewrite_replicas
! thermodynamical quantities
      real(kind=8),allocatable :: C_av(:),chi_M(:,:),chi_Q(:,:)
! errors on the different quantities
      real(kind=8),allocatable :: E_err_av(:),M_err_av(:,:)
! sums
      real(kind=8),allocatable :: M_sq_sum_av(:,:),E_sum_av(:),E_sq_sum_av(:),Q_sq_sum_av(:),Qp_sq_sum_av(:),Qm_sq_sum_av(:)
! energy and so on
      real(kind=8),allocatable :: E_av(:),qeulerp_av(:),qeulerm_av(:),M_av(:,:)
      real(kind=8),allocatable :: M_sum_av(:,:),vortex_av(:,:),chi_l(:,:)
      real(kind=8) :: spin_sum(4,shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5))
! table that stores the replicas
      real(kind=8), allocatable :: replicas(:,:,:,:,:,:)
! relaxation informations
      real(kind=8), allocatable :: relax(:,:,:)
! dummy slopes
      integer :: i,j,k,istart,istop
      integer :: ierr


#ifdef CPP_MPI
      include 'mpif.h'
#endif

! initialized the size of the tables
      size_table=n_Tsteps

#ifdef CPP_MPI
      size_M_replica=nRepProc
      if (irank_working.eq.0) then
       write(6,'(/,a,I6,a,/)') "setup from parallel tempering"
       write(6,'(a,I6,a,/)') "you are calculating",size_table," temperatures"
      endif
#else
      size_M_replica=n_Tsteps
      write(6,'(/,a,I6,a,/)') "setup from parallel tempering"
      write(6,'(/,a,I6,a,/)') "you are calculating",size_table," temperatures"
#endif

! allocate the temperatures
      allocate(kT_all(size_table),kT_updated(size_table))
      allocate(image_temp(size_table))
! allocate the matrix for up and down images
      allocate(nup(size_table),ndown(size_table),label_world(size_table),success_PT(size_table))
      allocate(E_av(size_table))
!     Errors for energy and magnetization
      allocate(E_err_av(size_table),M_err_av(3,size_table))
!     Energysum and Magnetizationsum and their squaresum
!     specific heat and suszeptibility
      allocate(C_av(size_table),chi_Q(4,size_table),Q_sq_sum_av(size_table),Qp_sq_sum_av(size_table),Qm_sq_sum_av(size_table))
!     magnetisation
      allocate(M_sum_av(3,size_table),M_sq_sum_av(3,size_table),chi_l(3,size_table),chi_M(3,size_table),M_av(3,size_table))
      allocate(vortex_av(3,size_table),qeulerp_av(size_table),qeulerm_av(size_table))
      allocate(E_sum_av(size_table),E_sq_sum_av(size_table),E_temp(size_table))

      allocate(relax(18,n_sizerelax,size_M_replica))
! allocate the matrix of the replicas
      allocate(replicas(shape_spin(1),shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5),size_M_replica))
! initialize the matrices
      do i_image=1,size_M_replica
       do i_m=1,shape_spin(5)
        do i_z=1,shape_spin(4)
         do i_y=1,shape_spin(3)
          do i_x=1,shape_spin(2)
           replicas(:,i_x,i_y,i_z,i_m,i_image)=spin(:,i_x,i_y,i_z,i_m)
          enddo
         enddo
        enddo
       enddo
      enddo

! open the data file where the thermodynamical quantities should be written
#ifdef CPP_MPI
    if (irank_working.eq.0) then
#endif
      OPEN(7,FILE='EM.dat',action='write',status='unknown', &
        & position='append',form='formatted')
        Write(7,'(29(a,15x))') '# 1:T','2:E_av','3:E_err','4:C','5:M','6:Mx','7:My','8:Mz', &
        & '9:M_err_x','10:M_err_y','11:M_err_z','12:chi_x','13:chi_y','14:chi_z','15:vx', &
        & '16:vy','17:vz','18:qeuler','19:ChiQ','20:Q+','21:ChiQ+','22:Q-','23:ChiQ-','24:Im','25:Qm','26:l_x','27:l_y','28:l_z','29:ChiQpQm'
#ifdef CPP_MPI
    endif
#endif

! open the data where the parallel tempering quantities will be written
#ifdef CPP_MPI
      if (irank_working.eq.0) then
#endif
        open(21,FILE='T-range.dat',action='write',status='unknown', &
        & position='append',form='formatted')

        if (.not.i_optTset) then
          write(21,'(a)') '#  1:T  2:fup(T)  3:A(T)'
        else
          write(21,'(a)') '#  1:Told   2:Tnew   3:fup(T)    4:D(T)    5:dfup(T)    6:A(T)    7:eta(T)'
        endif

#ifdef CPP_MPI
      endif
#endif

! Measured data
      E_av=0.0d0
      C_av=0.0d0
      M_av=0.0d0
      chi_Q=0.0d0
      E_err_av=0.0d0
      E_sum_av=0.0d0
      M_err_av=0.0d0
      M_sum_av=0.0d0
      vortex_av=0.0d0
      qeulerp_av=0.0d0
      qeulerm_av=0.0d0
      M_sq_sum_av=0.0d0
      E_sq_sum_av=0.0d0
      Q_sq_sum_av=0.0d0
      Qp_sq_sum_av=0.0d0
      Qm_sq_sum_av=0.0d0
      chi_l=0.0d0
      chi_M=0.0d0
      spin_sum=0.0d0

! temperature
      kT_all=0.0d0
      kT=0.0d0
      image_temp=(/(i,i=1,size_table)/) ! at the beginning, the first image has the first temperature, the second image, the second ...

! initialization of the variables
! convergence of the MC
      rate=0.0d0
      cone=coni
      Nsuccess=0.0d0
      tries=0.0d0

! control of the parallel tempering
      label_world=0.0d0
      nup=0.0d0
      ndown=0.0d0
      success_PT=0.0d0

! local variables
      qeulerp=0.0d0
      qeulerm=0.0d0
      vortex=0.0d0
      magnetization=0.0d0
      E_total=0.0d0
      E_decompose=0.0d0
      autocor_steps=T_auto
      kT=0.0d0
      E_temp=0.0d0

! relaxation variables
      relax=0.0d0
      relaxation_steps=n_sizerelax
      total_relax_steps=relaxation_steps*(2**N_temp-1)
      freq_rw=total_relax_steps/n_sizerelax

      istart=0
      istop=0
#ifdef CPP_MPI
! transfer data in the case of MPI
      istart=1+irank_working*nRepProc
      istop=(irank_working+1)*nRepProc
      image_temp(nRepProc+1:)=0
#endif

! restart of the parallel tempering
      i_reread_paratemp=.False.
      i_rewrite_paratemp=.False.
      i_reread_data=.False.
      i_rewrite_data=.False.
      i_reread_replicas=.False.
      i_rewrite_replicas=.False.

#ifdef CPP_MPI
    if (irank_working.eq.0) then
#endif
      write(6,'(/,a)') 'parallel tempering selected'
      write(6,'(a,I10,a,I10,/)') 'the number of relaxation steps will go from', relaxation_steps, ' to ', relaxation_steps*2**(N_temp-1)
#ifdef CPP_MPI
    endif
#endif

! initialize temperatures
#if CPP_MPI
     call ini_temp(kt_all,kTfin,kTini,size_table,irank_working,nRepProc,i_print_W)
#else
     call ini_temp(kt_all,kTfin,kTini,size_table,i_print_W)
#endif

! in case of restart of the temperature set
     call check_restart_read('temperature.in',kt_all,size_table)

     kT_updated=kt_all

#ifdef CPP_MPI

! kt_all on rank 0 must contain all the temperatures in the good order
     !call mpi_gather(kt_all,nRepProc,MPI_REAL8,kt_all(1:nRepProc),nRepProc,MPI_REAL8,0,MPI_COMM,ierr)
     kt_all=gather(kt_all(1:nRepProc),nRepProc,size_table,0,MPI_COMM)
#endif

! in case reinitialization of the fraction or writting out the fraction
     call check_restart_read('fraction.in',i_reread_paratemp)
     call check_restart_read('fraction.out',i_rewrite_paratemp)
! in case the measured data should be restarted from file
     call check_restart_read('measured_data.in',i_reread_data)
     call check_restart_read('measured_data.out',i_rewrite_data)
! in case the replicas should be restarted from file
     call check_restart_read('replicas.in',i_reread_replicas)
     call check_restart_read('replicas.out',i_rewrite_replicas)


! read the replicas from file in case it is necessary
#ifdef CPP_MPI
#else
     if (i_reread_replicas) call check_restart_read('replicas.in',replicas,shape_spin,size_table)
#endif

     do j_optset=1,N_temp

! reset the variables to 0 when the temperature set is changed
! Measured data
       E_av=0.0d0
       C_av=0.0d0
       chi_Q=0.0d0
       E_err_av=0.0d0
       E_sum_av=0.0d0
       M_err_av=0.0d0
       M_sum_av=0.0d0
       vortex_av=0.0d0
       qeulerp_av=0.0d0
       qeulerm_av=0.0d0
       M_sq_sum_av=0.0d0
       E_sq_sum_av=0.0d0
       Q_sq_sum_av=0.0d0
       chi_l=0.0d0
       chi_M=0.0d0
       nup=0.0d0
       ndown=0.0d0
       Nsuccess=0.0d0
       label_world=0.0d0
       spin_sum=0.0d0

       if (i_reread_paratemp) then
           call check_restart_read('fraction.in',nup,ndown,success_PT,label_world,image_temp,size_table)
           i_reread_paratemp=.False.
       endif

       if (i_reread_data) then
           call check_restart_read('measured_data.in',E_sum_av,E_sq_sum_av,E_err_av, &
      &   M_sum_av,M_sq_sum_av,M_err_av, qeulerp_av,qeulerm_av,Q_sq_sum_av,vortex_av, &
      &   size_table)
           i_reread_paratemp=.False.
       endif

!         T_relax is probably larger then T_relax_temp, this is because
!         one might want to thermalize the first iteration more as compared with the j_optset>2

       if (j_optset.gt.1) autocor_steps=T_relax_temp

! Relaxation of the System - SWAPPING!!
       do i_relax=1,relaxation_steps
! do loop over the temperatures
         do i_temp=1,size_M_replica

! find the replica which has this temperature
            i_image=image_temp(i_temp)
! load the actual temperature into kT
            kt=kt_updated(i_temp)
! initializing the variables above
            call DeriveValue(N_cell,replicas(:,:,:,:,:,i_image),shape_spin,tableNN,shape_tableNN,masque,shape_masque,indexNN,shape_index,E_decompose, &
      &   i_four,i_dm,i_biq,i_dip,i_stone,EA,h_ext)
            E_total=sum(E_decompose)

            Call CalculateAverages(replicas(:,:,:,:,:,i_image),shape_spin,masque,shape_masque,qeulerp,qeulerm,vortex,magnetization,n_system)

            Do i_MC=1,autocor_steps*N_cell
                Call MCStep(state,E_total,E_decompose,Magnetization,kt,Nsuccess,rate,tries,cone,n_system, &
            & replicas(:,:,:,:,:,i_image),shape_spin,tableNN,shape_tableNN,masque,shape_masque,indexNN,shape_index,h_ext,EA, &
            & i_biq,i_dip,i_DM,i_four,i_stone,ising,i_print_W,equi,overrel,sphere,underrel,world)
            enddo

! In case T_relax set to zero at least one MCstep is done
            Call MCStep(state,E_total,E_decompose,Magnetization,kt,Nsuccess,rate,tries,cone,n_system, &
       & replicas(:,:,:,:,:,i_image),shape_spin,tableNN,shape_tableNN,masque,shape_masque,indexNN,shape_index,h_ext,EA, &
       & i_biq,i_dip,i_DM,i_four,i_stone,ising,i_print_W,equi,overrel,sphere,underrel,world)

! calculate the topocharge
            call topo(replicas(:,:,:,:,:,i_image),shape_spin,masque,shape_masque,qeulerp,qeulerm)

! CalculateAverages makes the averages from the sums
            Call CalculateAverages(qeulerp_av(i_temp),qeulerm_av(i_temp),Q_sq_sum_av(i_temp),Qp_sq_sum_av(i_temp),Qm_sq_sum_av(i_temp) &
       &  ,vortex_av(:,i_temp),vortex &
       &  ,E_sum_av(i_temp),E_sq_sum_av(i_temp),M_sum_av(:,i_temp),M_sq_sum_av(:,i_temp),E_total,Magnetization,spin_sum &
       &  ,replicas(:,:,:,:,:,i_image),shape_spin,masque,shape_masque)

!           save the data for thermalization

            if (print_relax) then

               i_pos=i_relax+n_sizerelax*(2**(j_optset-1)-1)

               pos=dble(i_pos)
               i_rw=mod(i_pos,freq_rw)

               if (i_rw.eq.0) call store_relaxation(Relax,i_temp,n_sizerelax,i_pos/freq_rw,pos, &
       &  sum(E_decompose)/dble(N_cell),E_decompose,dble(N_cell),kt,Magnetization,rate,cone,qeulerp,qeulerm)

            endif

! save the total energy of each replicas
         E_temp(i_temp)=E_total

         enddo     ! enddo of the images


#ifdef CPP_MPI
         if (i_ghost) then
!             call paratemp(E_total,kT,state,label &
!           & ,nup,ndown,i_thousand,Nsuccess,n_thousand, &
!           & vortex_mpi,qeulerp_mpi,qeulerm_mpi,kt_mpi,qeulerp,qeulerm,E_mpi,E_sq_mpi,M_sq_mpi, &
!           & Magnetization,M,vortex,j_optset,N_temp,n_ghost)
         else
! give the knowledge of all the temperature and all the energies to all the processors
            call paratemp_gather(istart,istop,irank_working,nRepProc, &
          &   size_table,MPI_COMM,kt_updated,E_temp,image_temp, &
          &   qeulerp_av,qeulerm_av,Q_sq_sum_av,Qp_sq_sum_av,Qm_sq_sum_av, &
          &   vortex_av,E_sum_av,E_sq_sum_av,M_sum_av,M_sq_sum_av, &
          &   label_world)

         endif
#endif

#ifdef CPP_MPI
         if (irank_working.eq.0) then
#endif
         call paratemp(state,label_world,nup,ndown,i_relax,success_PT,relaxation_steps, &
      &    kt_updated,image_temp,E_temp,N_temp,size_table,kt_all)
#ifdef CPP_MPI
         endif

! reorder the temperature and the replicas
! The replicas are actually attached to the processors so they stay remote and can not be rearranged
! this part makes sure that the temperature are on the processors of the replicas and that Image_temp points to the right temperature
! and then scatter the data on the good remote processors
         call paratemp_scatter(isize,irank_working,nRepProc,kt_updated,image_temp,MPI_COMM, &
          &   qeulerp_av,qeulerm_av,Q_sq_sum_av,Qp_sq_sum_av,Qm_sq_sum_av, &
          &   vortex_av,E_sum_av,E_sq_sum_av,M_sum_av,M_sq_sum_av, &
          &   label_world)

#endif
         E_temp=0.0d0
      enddo  !  end of the relaxation steps (swap-steps)


    if (gra_log) then
     !write the Spinsefiles
        do i_temp=1,size_M_replica
! load the actual temperature into kT
           call CreateSpinFile(kt_updated(i_temp)/k_b,replicas(:,:,:,:,:,i_temp),shape_spin)
        enddo
    endif


!        write the Equilibrium-Files in case of paratemp

#ifdef CPP_MPI
      if (print_relax) call write_relaxation(Relax,kT_updated,irank_working,n_sizerelax,size_M_replica,size_table,Cor_log)
#else
      if (print_relax) call write_relaxation(Relax,kT_all,n_sizerelax,size_table,Cor_log)
#endif

! at the end of the relaxation and the image loops, calculate the thermodynamical quantitites
#ifdef CPP_MPI
       call paratemp_gather(istart,istop,irank_working,nRepProc, &
          &   size_table,MPI_COMM,kt_updated,E_temp,image_temp, &
          &   qeulerp_av,qeulerm_av,Q_sq_sum_av,Qp_sq_sum_av,Qm_sq_sum_av, &
          &   vortex_av,E_sum_av,E_sq_sum_av,M_sum_av,M_sq_sum_av, &
          &   label_world)

       if (irank_working.eq.0) then
#endif
       call Calculate_thermo(Cor_log,relaxation_steps,dble(N_cell),kT_updated,E_sq_sum_av,E_sum_av,M_sq_sum_av, &
    &    C_av,chi_M,E_av,E_err_av,M_err_av,qeulerp_av,qeulerm_av,vortex_av,Q_sq_sum_av,Qp_sq_sum_av,Qm_sq_sum_av,chi_Q,chi_l, &
    &    M_sum_av,M_av,size_table)

! write the thermodynamical quantities in EM.dat

        do i=1,size_table
            Write(7,'(24(E20.10E3,2x),3(E20.10E3))') kT_updated(i)/k_B ,E_av(i), E_err_av(i), C_av(i), norm(M_sum_av(:,i))/N_cell/relaxation_steps, &
   &     M_av(:,i), M_err_av(:,i), chi_M(:,i), vortex_av(:,i), qeulerm_av(i)+qeulerp_av(i), chi_Q(1,i), &
   &     qeulerp_av(i), chi_Q(2,i), qeulerm_av(i), chi_Q(3,i), chi_l(:,i), chi_Q(4,i)
        enddo
         Write(7,'(/)')

! calculate the current, diffusivity... for the paratemp
      call calculate_diffusion(kT_all,kt_updated,nup,ndown,success_PT,relaxation_steps,size_table,i_optTset)
#ifdef CPP_MPI
      endif
#endif
      if ((i_optTset).and.(j_optset.ne.N_temp)) then

! update the temperature set
      kT_all=kt_updated

      write(6,'(/,a,I5,a,I5)') '#-------- parallel tempering step', j_optset ,' of ', N_temp
      write(6,'(a,I10,a,I10,/)') '#-------- number of relaxation steps goes from', relaxation_steps ,' --- to --- ', relaxation_steps*2

      relaxation_steps=relaxation_steps*2

#ifdef CPP_MPI
      if (irank_working.eq.0) write(6,'(/,a,/)') 'reset of the current and the labels'
#else
      write(6,'(/,a,/)') 'reset of the current and the labels'
#endif

      endif

! the data for paratemp in case a restart is needed
    if (i_rewrite_paratemp) call check_restart_write('fraction.out',nup,ndown,success_PT,label_world,image_temp,size_table)
    if (i_rewrite_data) call check_restart_write('measured_data.out',E_sum_av,E_sq_sum_av,E_err_av, &
      &   M_sum_av,M_sq_sum_av,M_err_av, qeulerp_av,qeulerm_av,Q_sq_sum_av,vortex_av, &
      &   size_table)
! write the replicas from file in case it is necessary
     if (i_rewrite_replicas) call check_restart_write('replicas.out',replicas,shape_spin,size_table)

    enddo !over j_MC->N_tempcd


!    get ready to start the next calculation with the new temperature set
!    exists = .False.
!    inquire (file='restart',exist=exists) !attention: Be carefull with reading old Spinstructures at a restart!
#ifdef CPP_MPI
    if (irank_working.eq.0) call check_restart_write('temperature.out',kt_updated,size_table)
#else
       call check_restart_write('temperature.out',kt_updated,size_table)
#endif

#ifdef CPP_MPI
        if (irank_working.eq.0) write(6,'(a)') 'parallel tempering is done'
#else
        write(6,'(a)') 'parallel tempering is done'
#endif

!close all the file
    close(21)
    close(7)

      end subroutine parallel_tempering
