!
! Routine that does the Monte Carlo (and not the parallel tempering)
!
      subroutine montecarlo(N_cell,n_system,state,world,n_relaxation,n_sizerelax, &
            &    spin,shape_spin,tableNN,shape_tableNN,masque,shape_masque,indexNN,shape_index, &
            &    i_qorien,i_biq,i_dip,i_DM,i_four,i_stone,ising,i_print_W,equi,overrel,sphere,underrel,i_ghost, &
            &    Periodic_log,gra_topo,CalEnergy,CalTheta,Gra_log,spstmL,gra_fft,&
            &    i_separate,i_average,i_topohall,print_relax,Cor_log, &
            &    n_Tsteps,coni,Total_MC_Steps,T_auto,EA,T_relax,kTfin,kTini,h_ext, &
            &    n_ghost,nRepProc)     ! MPI variable
      use m_constants, only : k_b,pi
      use m_vector, only : norm
      use m_Corre
      use mtprng
      use m_check_restart
      use m_topocharge
      use m_qorien
      use m_check_restart
      use m_average_MC
      use m_write_spin
      use m_set_temp
      use m_topocharge_all
      use m_mpi_prop, only : isize,irank_working
      use m_createspinfile
#ifdef CPP_MPI
      use m_mpi_prop, only : MPI_COMM
      use m_gather_reduce
#endif
      implicit none
! interface for the relaxation and the MCsteps routine
      interface
        subroutine Relaxation(N_cell,n_system,kT,state,E_total,E_decompose,Magnetization,qeulerp,qeulerm,vortex, &
      &    n_relaxation,n_sizerelax,T_relax,acc,rate,tries,cone,print_relax,h_ext,EA, &
      &    spin,shape_spin,tableNN,shape_tableNN,masque,shape_masque,indexNN,shape_index, &
      &    i_biq,i_dip,i_DM,i_four,i_stone,ising,i_print_W,equi,overrel,sphere,underrel,n_world,i_ghost)
          use mtprng
          integer, intent(in) :: shape_index(2),shape_spin(5),shape_tableNN(6),shape_masque(4)
          real(kind=8), intent(inout) :: spin(shape_spin(1),shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5))
          type(mtprng_state), intent(inout) :: state
          real(kind=8), intent(inout) :: qeulerp,qeulerm,vortex(3),cone,acc,rate,tries
          real(kind=8), intent(inout) :: E_total,magnetization(3),E_decompose(8)
          real(kind=8), intent(in) :: kT,h_ext(3),EA(3)
          integer, intent(in) :: tableNN(shape_tableNN(1),shape_tableNN(2),shape_tableNN(3),shape_tableNN(4),shape_tableNN(5),shape_tableNN(6))
          integer, intent(in) :: masque(shape_masque(1),shape_masque(2),shape_masque(3),shape_masque(4))
          integer, intent(in) :: indexNN(shape_index(1),shape_index(2)),n_world
          integer, intent(in) :: n_relaxation,T_relax,N_cell,n_sizerelax,n_system
          logical, intent(in) :: print_relax,i_biq,i_dip,i_DM,i_four,i_stone,ising,i_print_W,equi,overrel,sphere,underrel,i_ghost
        end subroutine

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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!variable part
      integer, intent(in) :: N_cell,n_system,n_relaxation,n_sizerelax,n_Tsteps,Total_MC_Steps,T_auto,T_relax,n_ghost,nRepProc
      integer, intent(in) :: shape_spin(5),shape_tableNN(6),shape_masque(4),shape_index(2),world
      integer, intent(in) :: tableNN(shape_tableNN(1),shape_tableNN(2),shape_tableNN(3),shape_tableNN(4),shape_tableNN(5),shape_tableNN(6))
      integer, intent(in) :: masque(shape_masque(1),shape_masque(2),shape_masque(3),shape_masque(4))
      integer, intent(in) :: indexNN(shape_index(1),shape_index(2))
      logical, intent(in) :: i_qorien,Periodic_log(3),gra_topo,CalEnergy,CalTheta,Gra_log,spstmL,gra_fft,i_separate,i_average, &
    &   i_topohall,print_relax,Cor_log,i_biq,i_dip,i_DM,i_four,i_stone,ising,i_print_W,equi,overrel,sphere,underrel,i_ghost
      real(kind=8), intent(in) :: kTfin,kTini,coni,EA(3),h_ext(3)
      real(kind=8), intent(inout) :: spin(shape_spin(1),shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5))
      type(mtprng_state), intent(inout) :: state
! internal variables
! slope of the MC
      integer :: i_relax,n_kT,n_MC
!restart variables
      integer :: restart_MC_steps
      logical :: i_restart
! size of the table for the mpi averaging
      integer :: size_table
! variable for the temperature
      real(kind=8) :: kT
! variables that being followed during the simulation
      real(kind=8) :: qeulerp,qeulerm,vortex(3),magnetization(3),E_total
! contribution of the different energy parts
      real(kind=8) :: E_decompose(8)
! thermodynamical quantities
      real(kind=8),allocatable :: C_av(:),chi_M(:,:),chi_Q(:,:)
! errors on the different quantities
      real(kind=8),allocatable :: E_err_av(:),M_err_av(:,:)
! sums
      real(kind=8),allocatable :: M_sq_sum_av(:,:),E_sum_av(:),E_sq_sum_av(:),Q_sq_sum_av(:),Qp_sq_sum_av(:),Qm_sq_sum_av(:)
! energy and so on
      real(kind=8),allocatable :: E_av(:),qeulerp_av(:),qeulerm_av(:),kt_all(:),M_av(:,:)
      real(kind=8),allocatable :: M_sum_av(:,:),vortex_av(:,:),chi_l(:,:)
      real(kind=8) :: spin_sum(4,shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5))
      real(kind=8) :: angle_sum(2,shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5))

! variable for the convergence of the MC
      real(kind=8) :: acc,rate,tries,cone
      character(len=30) :: fname,toto,fname2
      integer :: i,i_x,i_y,i_z,i_m,k

! dummys
      integer :: ierr
#ifdef CPP_MPI

      if (i_separate) size_table=isize*nRepProc
!      if (i_ghost) size_table=isize/n_ghost
      if (irank_working.eq.0) write(6,'(/,a,I6,a,/)') "you are calculating",size_table," temperatures"
#else
      size_table=n_Tsteps
      write(6,'(/,a,I6,a,/)') "you are calculating",size_table," temperatures"
#endif

      allocate(E_av(size_table))
!     Errors for energy and magnetization
      allocate(E_err_av(size_table),M_err_av(3,size_table))
!     Energysum and Magnetizationsum and their squaresum
!     specific heat and suszeptibility
      allocate(C_av(size_table),chi_Q(4,size_table))
! everything for the topological charge
      allocate(Q_sq_sum_av(size_table),Qp_sq_sum_av(size_table),Qm_sq_sum_av(size_table))
!     magnetisation
      allocate(M_sum_av(3,size_table),M_sq_sum_av(3,size_table),chi_l(3,size_table),chi_M(3,size_table),M_av(3,size_table))
      allocate(vortex_av(3,size_table),qeulerp_av(size_table),qeulerm_av(size_table))
      allocate(kt_all(size_table),E_sum_av(size_table),E_sq_sum_av(size_table))

! updating data during the simulation
      qeulerp=0.0d0
      qeulerm=0.0d0
      vortex=0.0d0
      magnetization=0.0d0
      E_total=0.0d0
      E_decompose=0.0d0
      kt_all=0.0d0

! initializing the variables above
      call DeriveValue(N_cell,spin,shape_spin,tableNN,shape_tableNN,masque,shape_masque,indexNN,shape_index,E_decompose, &
           &   i_four,i_dm,i_biq,i_dip,i_stone,EA,h_ext)
      E_total=sum(E_decompose)

      Call CalculateAverages(spin,shape_spin,masque,shape_masque,qeulerp,qeulerm,vortex,magnetization,n_system)

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
      Qp_sq_sum_av=0.0d0
      Qm_sq_sum_av=0.0d0
      chi_l=0.0d0
      chi_M=0.0d0
      M_av=0.0d0

! statistics on the MC
      acc=0.0d0
      rate=0.0d0
      tries=0.0d0
      cone=coni

! initialization of the spin_sum and angle_sum. Used in case of the FFT
      spin_sum=0.0d0
      angle_sum=0.0d0

! values in case of restart
      restart_MC_steps=0
      i_restart=.False.

! initialize the temperatures
#ifdef CPP_MPI
      call ini_temp(kt_all,kTfin,kTini,size_table,irank_working,nRepProc,i_print_W)
#else
      call ini_temp(kt_all,kTfin,kTini,size_table,i_print_W)
#endif

#ifdef CPP_DEBUG
      write(*,*) irank_working,kt_all/k_B
#endif

!      inquire (file='restart',exist=i_restart)
!      if (i_restart) then
!         write(6,'(/,a,/)') 'restart from previous configurations'
!         call check_restart_read(spin,shape_spin,i_separate,i_paratemp,size_table,kt &
!     &     ,E_av,E_err_av,M_err_av,Q_sq_sum_av,M_sum_av,qeulerp_av,qeulerm_av,vortex_av,kt_all,E_sum_av,E_sq_sum_av)
!      endif


      Do n_kT=1,n_Tsteps

       kt=kt_all(n_kT)

        call Relaxation(N_cell,n_system,kT,state,E_total,E_decompose,Magnetization,qeulerp,qeulerm,vortex, &
      &    n_relaxation,n_sizerelax,T_relax,acc,rate,tries,cone,print_relax,h_ext,EA, &
      &    spin,shape_spin,tableNN,shape_tableNN,masque,shape_masque,indexNN,shape_index, &
      &    i_biq,i_dip,i_DM,i_four,i_stone,ising,i_print_W,equi,overrel,sphere,underrel,world,i_ghost)

!       Monte Carlo steps, calculate the values

            do n_MC=1+restart_MC_steps,Total_MC_Steps+restart_MC_steps

!         Monte Carlo steps for independency

                Do i_relax=1,T_auto*N_cell

                    Call MCstep(state,E_total,E_decompose,Magnetization,kt,acc,rate,tries,cone,n_system, &
            & spin,shape_spin,tableNN,shape_tableNN,masque,shape_masque,indexNN,shape_index,h_ext,EA, &
            & i_biq,i_dip,i_DM,i_four,i_stone,ising,i_print_W,equi,overrel,sphere,underrel,world)

                End do

                Call MCstep(state,E_total,E_decompose,Magnetization,kt,acc,rate,tries,cone,n_system, &
            & spin,shape_spin,tableNN,shape_tableNN,masque,shape_masque,indexNN,shape_index,h_ext,EA, &
            & i_biq,i_dip,i_DM,i_four,i_stone,ising,i_print_W,equi,overrel,sphere,underrel,world)

! Calculate the topological charge and the vorticity
                call topo(spin,shape_spin,masque,shape_masque,qeulerp,qeulerm)

! CalculateAverages makes the averages from the sums
                Call CalculateAverages(qeulerp_av(n_kT),qeulerm_av(n_kT),Q_sq_sum_av(n_kT),Qp_sq_sum_av(n_kT),Qm_sq_sum_av(n_kT),vortex_av(:,n_kT),vortex &
                &  ,E_sum_av(n_kT),E_sq_sum_av(n_kT),M_sum_av(:,n_kT),M_sq_sum_av(:,n_kT),E_total,Magnetization,spin_sum,spin,shape_spin, &
                    masque,shape_masque)

               if (Cor_log) chi_l(:,n_kT)=chi_l(:,n_kT)+Correlation(spin_sum,spin(4:6,:,:,:,:),shape_spin,n_MC,dble(N_cell))


! Calculate the sum of the spin components and angles for average
                if (CalTheta) Call SphericalCoordinates(spin(4:6,:,:,:,:),shape_spin,angle_sum)

!**************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!
            end do ! over n_MC
!***************************
!!!!!!!!!!!!!***************************!!!!!!
!!!!!!!!!!!!!***************************!!!!!!

! FFT calculation
!       if (total_MC_steps.ne.0) call MC_fft(sum(spin_sum(1:3,:,:,:,:),5)/dble(total_MC_steps), &
!     & spin(1:3,:,:,:,1),Im,Qim,dim_lat,N_cell)

!#ifdef CPP_MPI
!
!        if (i_ghost.or.i_separate.or.i_average) then
!            if (gra_topo) then
!                call mpi_average(n_kt,vortex_mpi,qeulerp_mpi, &
!                & qeulerm_mpi,map_vort,mapvort_mpi,mapeul_mpi,map_toto, &
!                & Im_mpi,Qim_mpi,Im,Qim,qeulerm,qeulerp,vortex, &
!                & E_sum,E_sq_sum,M_sum,M_sq_sum,C,chi,E_err,M,E_av,M_err,dble(N_cell))
!            else
!                if (total_MC_steps.ne.0) then
!                    if (i_ghost) then
!                        call mpi_average(vortex_mpi(:,1),qeulerp_mpi(1), &
!                        & qeulerm_mpi(1),Im_mpi(1,:),Qim_mpi(1,:),Im,Qim,qeulerm,qeulerp,vortex, &
!                        & E_sum,E_sq_sum,M_sum,M_sq_sum,C(1),chi(1),E_err(1),M(:,1),E_av(1),M_err(1),n_ghost,MPI_COMM_BOX)
!                    else
!                        call mpi_average(n_kt,vortex_mpi,qeulerp_mpi, &
!                        & qeulerm_mpi,Im_mpi,Qim_mpi,Im,Qim,qeulerm,qeulerp,vortex, &
!                        & E_sum,E_sq_sum,M_sum,M_sq_sum,C,chi,E_err,M,E_av,M_err,dble(N_cell))
!                    endif
!                endif
!            endif
!       endif
!
!#endif

       if (n_Tsteps.ne.0) call Calculate_thermo(Cor_log,total_MC_steps, &
    &    dble(N_cell),kT_all(n_kt),E_sq_sum_av(n_kt),E_sum_av(n_kt),M_sq_sum_av(:,n_kt), &
    &    C_av(n_kt),chi_M(:,n_kt),E_av(n_kt),E_err_av(n_kt),M_err_av(:,n_kt),qeulerp_av(n_kt),qeulerm_av(n_kt),vortex_av(:,n_kt), &
    &    Q_sq_sum_av(n_kt),Qp_sq_sum_av(n_kt),Qm_sq_sum_av(n_kt), &
    &    chi_Q(:,n_kt),chi_l(:,n_kt), &
    &    M_sum_av(:,n_kt),M_av(:,n_kt))


        write(6,'(5(a,f18.9,2x))') 'M= ',norm(M_av(:,n_kt)), &
     & 'E= ',E_av(n_kT),'Q+= ',qeulerp_av(n_kT),'Q-= ',qeulerm_av(n_kT),'Q= ',qeulerp_av(n_kT)+qeulerm_av(n_kT)


       if (Gra_log) then

            Call WriteSpinAndCorrFile(kt,spin,shape_spin)

            call CreateSpinFile(kt/k_B,spin,shape_spin)

       endif

!ccccccccccccccccccccccccccccccccccccc
! Calculate the topological charge
!cccccccccccccccccccccccccccccccccccc
        if (gra_topo) then
           if (world.eq.2) then
               if (shape_spin(5).eq.1) then
                  call topo_map(spin(4:6,:,:,1,1),shape_spin,kt/k_B,gra_topo,Periodic_log,i_separate)
               else
                  call topo_map(spin(4:6,:,:,1,:),shape_spin,kt/k_B,gra_topo,Periodic_log,i_separate)
               endif
           elseif (world.eq.1) then
              write(6,'(a)') 'topological graphs not coded for 1D'
           else
              write(6,'(a)') 'topological graphs not coded for 0D or 3D system'
           endif
        endif

!ccccccccccccccccccccccccccccccccccccc
! Calculate the oriented topological charge
!cccccccccccccccccccccccccccccccccccc
        if (i_qorien) then
           if (world.eq.2) then
               if (shape_spin(5).eq.1) then
                  call qorien(spin(4:6,:,:,1,1),shape_spin)
               else
                  call qorien(spin(4:6,:,:,1,:),shape_spin)
               endif
           elseif (world.eq.1) then
             write(6,'(a)') 'spacial resolution of the gauge not coded for 1D'
           else
             write(6,'(a)') 'topological graphs not coded for 0D or 3D system'
           endif
        endif

! ===========================================================
! Write the averaged coordinate of spin in spherical coordinates
! ===========================================================
99999 format(5f14.8)

        if (CalTheta) then

           Write(fname,'(f8.4)') kT/k_B
           toto = Trim(Adjustl(AdjustR(fname)))

            write(fname,'(a,14a,a)')'Theta_T_',(toto(i:i),i=1,len_trim(toto)),'.dat'

            open(666,file =Trim(Adjustl(fname)), form = 'formatted', &
            &    status ='unknown',  action = 'write')

            do i_m=1,shape_spin(5)
                do i_z=1,shape_spin(4)
                    do i_y=1,shape_spin(3)
                        do i_x=1,shape_spin(2)

                            Write(666,99999) (Spin(k,i_x,i_y,i_z,i_m),k=1,3), &
                            &    (angle_sum(k,i_x,i_y,i_z,i_m)/dble(total_MC_steps),k=1,2)
                        enddo
                    enddo
                enddo
            enddo

            close(666)

        endif
! ===========================================================
! Computation of the energy density on the lattice
! ===========================================================

        if (CalEnergy) then

            write(fname,'(a,14a,a)')'Energy_T_',(toto(i:i),i=1, &
            &  len_trim(toto)),'.dat'
            write(fname2,'(a,14a,a)')'DensityOfEnergy_T_',(toto(i:i),i=1, &
            &  len_trim(toto)),'.dat'

            Call EnergyDensity(fname,fname2,spin,shape_spin,tableNN,shape_tableNN,masque,shape_masque,indexNN,shape_index)
        endif

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Calcul of the fft
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#ifndef CPP_MPI
!        if (gra_fft) call fft(shape_spin(2:4),kT)
#endif

        Write(6,'(I6,a,I6,a,f8.4,a,/)')  n_kT, ' nd step out of ', n_Tsteps,' steps. T=', kT/k_B,' Kelvin'

        end do !over n_kT
!--------------------------------------------------------------
!!!!!!!!!!!!!***************************!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!***************************!!!!!!!!!!!!!!!!!!!!!!!

! The program of Tobias is used only at last iteration
#ifdef CPP_MPI
        if ((irank_working.eq.0).AND.spstmL) call spstm
#else
        if (spstmL) call spstm
#endif

#ifdef CPP_MPI
    if (irank_working.eq.0) then
#endif
      OPEN(7,FILE='EM.dat',action='write',status='unknown', &
        & position='append',form='formatted')
        Write(7,'(27(a,15x))') '# 1:T','2:E_av','3:E_err','4:C','5:M','6:Mx','7:My','8:Mz', &
        & '9:M_err_x','10:M_err_y','11:M_err_z','12:chi_x','13:chi_y','14:chi_z','15:vx', &
        & '16:vy','17:vz','18:qeuler','19:Chi_q','20:Q+','21:Chi_qp','22:Q-','23:Chi_qm', &
        & '24:l_x','25:l_y','26:l_z','27:Chi_QpQm'
#ifdef CPP_MPI
    endif
#endif

#ifdef CPP_MPI

      if (i_separate) call end_gather(kT_all,E_av,E_err_av,C_av,M_sum_av,M_err_av,chi_M,vortex_av,chi_Q,qeulerp_av,qeulerm_av,chi_l, &
                          & size_table,irank_working,n_Tsteps,MPI_COMM)

   ! case of the ghost
!        if (i_ghost) then
!            n_Tsteps=isize/n_ghost
!            if (irank_box.eq.0) call end_gather(shape_masque,vortex_mpi,qeulerp_mpi,qeulerm_mpi,Im_mpi,Qim_mpi,kt_mpi,mpi_l, &
!            &     E_av(:),C(:),M_err(:),E_err(:),chi(:),M(:,:),n_Tsteps,MPI_COMM_MASTER)
!        endif

!        if ((i_separate.or.i_paratemp).and.(.not.i_ghost)) then
!            n_Tsteps=isize

!            call end_gather(shape_masque,vortex_mpi,qeulerp_mpi,qeulerm_mpi,Im_mpi,Qim_mpi,kt_mpi,mpi_l, &
!            &    E_av(:),C(:),M_err(:),E_err(:),chi(:),M(:,:),isize,MPI_COMM)
!        endif

!        if ((irank_working.eq.0).and.(.not.i_paratemp)) then
!            do i=1,n_Tsteps
!                Write(7,'(23(E20.10E3,2x),E20.10E3)') kt_mpi(i)/k_B,E_av(i),E_err(i),C(i), &
!                & norm(M(:,i)), M(:,i),M_err(:,i), chi(i), vortex_mpi(:,i) &
!                & ,qeulerp_mpi(i)+qeulerm_mpi(i),qeulerp_mpi(i),qeulerm_mpi(i), &
!                & Im_mpi(i,1),Qim_mpi(i,1),mpi_l(:,i)
!            enddo
!        endif
#endif

#ifdef CPP_MPI
      if (irank_working.eq.0) then

        do i=1,size_table
            Write(7,'(27(E20.10E3,2x),E20.10E3)') kT_all(i)/k_B ,E_av(i), E_err_av(i), C_av(i), norm(M_sum_av(:,i))/N_cell/(Total_MC_Steps+restart_MC_steps), &
     &             M_sum_av(:,i), M_err_av(:,i), chi_M(:,i), vortex_av(:,i), qeulerm_av(i)+qeulerp_av(i), chi_Q(1,i), &
     &             qeulerp_av(i), chi_Q(2,i), qeulerm_av(i), chi_Q(3,i), chi_l(:,i), chi_Q(4,i)
        enddo

        close(7) !Close EM.dat

      endif
#else
! write the data into a file
        do i=1,n_Tsteps
            Write(7,'(27(E20.10E3,2x),E20.10E3)') kT_all(i)/k_B ,E_av(i), E_err_av(i), C_av(i), norm(M_sum_av(:,i))/N_cell/(Total_MC_Steps+restart_MC_steps), &
     &             M_sum_av(:,i), M_err_av(:,i), chi_M(:,i), vortex_av(:,i), qeulerm_av(i)+qeulerp_av(i), chi_Q(1,i), &
     &             qeulerp_av(i), chi_Q(2,i), qeulerm_av(i), chi_Q(3,i), chi_l(:,i), chi_Q(4,i)
        enddo

      close(7) !Close EM.dat
#endif

    end subroutine montecarlo
