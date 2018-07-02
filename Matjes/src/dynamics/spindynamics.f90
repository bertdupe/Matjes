      subroutine spindynamics(state,spin,shape_spin,tableNN,shape_tableNN,masque,shape_masque,indexNN,shape_index,world, &
    &   N_cell, &
    &   i_biq,i_dm,i_four,i_dip,gra_topo,gra_log,gra_freq, &
    &   ktini,ktfin,EA,h_ext)
      use m_measure_temp
      use m_fieldeff
      use m_solver
      use m_dynamic
      use m_vector, only : cross,norm,norm_cross
      use m_sd_averages
      use m_randist
      use m_constants, only : pi,k_b,hbar
      use m_topo_sd
      use m_eval_Beff
      use mtprng
      use m_error_correction_SD
      use m_write_spin
      use m_energyfield
      use m_createspinfile
      use m_energy
      use m_local_energy
      use m_Torana
#ifndef CPP_BRUTDIP
      use m_setup_dipole, only : mmatrix
#endif
#ifdef CPP_MPI
      use m_parameters, only : i_ghost
      use m_mpi_prop, only : MPI_COMM,irank,isize,start
      use m_reconstruct_mat
#endif
      implicit none
! input
      integer, intent(in) :: shape_index(2),shape_spin(5),shape_tableNN(6),shape_masque(4),N_cell,world(:),gra_freq
      real(kind=8), intent(inout) :: spin(shape_spin(1),shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5))
      type(mtprng_state), intent(inout) :: state
      integer, intent(in) :: tableNN(shape_tableNN(1),shape_tableNN(2),shape_tableNN(3),shape_tableNN(4),shape_tableNN(5),shape_tableNN(6))
      integer, intent(in) :: masque(shape_masque(1),shape_masque(2),shape_masque(3),shape_masque(4)),indexNN(shape_index(1),shape_index(2))
      logical, intent(in) :: i_biq,i_dm,i_four,i_dip,gra_topo,gra_log
      real(kind=8), intent(in) :: ktini,ktfin,EA(3),h_ext(3)
! internal
      integer :: i,j,l,k,h,i_dum
      real(kind=8) :: spinafter(4,shape_tableNN(3),shape_tableNN(4),shape_tableNN(5),shape_tableNN(6))
      real(kind=8) :: Bini(3,shape_tableNN(3),shape_tableNN(4),shape_tableNN(5),shape_tableNN(6))
      real(kind=8) :: spinini(4,shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5))
      real(kind=8) :: dum_norm,qeuler,vortex(3),Mdy(3),Edy,stmtorquebp,Beff(3),check1,check2,Eold,check3
      real(kind=8) ::  qx,qy,qz,Mx,My,Mz,vx,vy,vz,check(2),test_torque,Einitial,ave_torque
      real(kind=8) :: dumy(4),ds(3),security(2),B(3),step(3),steptor(3),stepadia(3),stepsttor(3),steptemp(3)
      real(kind=8) :: timestep_ini,real_time,h_int(3)
      character(len=30) :: fname,toto
      integer :: i_x,i_y,i_z,i_m,iomp(3),N_site,ierr,shape_spinini(5)
! temperature (incase it is needed)
      real(kind=8) :: kT
! the computation time
      real(kind=8) :: computation_time
! parameter for the Heun integration scheme
      integer :: total_neigh,n_heun
      real(kind=8) :: maxh
! parameter for the rkky integration
      real(kind=8) :: path
      integer :: N_site_comm
! dumy
      logical :: said_it_once,i_anatorque
! starting and ending points of the sums
      integer :: Mstop
#ifndef CPP_MPI
      integer, dimension(3), parameter :: start=0
#endif
      integer :: Xstart,Xstop,Ystart,Ystop,Zstart,Zstop
#ifdef CPP_MPI
      real(kind=8) :: mpi_check(2),trans(3)

      include 'mpif.h'

      trans=0.0d0
#endif
! VERY IMPORTANT PART THAT DEFINES THE BOUNDARIES OF THE SUM
! starting point and ending points in the sums
      Xstart=start(1)+1
      Xstop=start(1)+shape_tableNN(3)
      Ystart=start(2)+1
      Ystop=start(2)+shape_tableNN(4)
      Zstart=start(3)+1
      Zstop=start(3)+shape_tableNN(5)
      Mstop=shape_tableNN(6)
      shape_spinini=shape(spinini)

      N_site_comm=(Xstop-Xstart+1)*(Ystop-Ystart+1)*(Zstop-Zstart+1)*shape_spin(5)
      i_anatorque=.False.

#ifdef CPP_MPI
      if (irank.eq.0) then
#endif
      OPEN(7,FILE='EM.dat',action='write',status='replace',form='formatted')
      Write(7,'(19(a,2x))') '# 1:time','2:E_av','3:M', &
     &  '4:Mx','5:My','6:Mz','7:vorticity','8:vx', &
     &  '9:vy','10:vz','11:qeuler','12:STMtorque','13:torque','14:T=', &
     &  '15:Tfin=','16:Ek=','17:Hx','18:Hy=','19:Hz='

! check the convergence
      open(8,FILE='convergence.dat',action='write',form='formatted')
#ifdef CPP_MPI
      endif
#endif
      call rw_dyna(shape_spin(2:4))
      timestep_ini=timestep


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Analysis of the Torques
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      inquire (file='analyse-Torque.in',exist=i_anatorque)
      if (i_anatorque) then
         write(6,'(/,a,/)') 'Analysis of the STT routine selected'
         call Torana(spin,shape_spin,indexNN,shape_index,masque,shape_masque,tableNN,shape_tableNN, &
              torque_AFL,torque_FL,damping,Ipol,EA,h_ext, &
              i_DM,i_four,i_biq,i_dip,i_torque,stmtorque, &
              timestep,adia,nonadia,storque,1.0d0)
      endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! End of analysis of the Torques
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      select case (integtype)
       case(3)
!-----------------------
! Heun
#ifdef CPP_MPI
        if (irank.eq.0) write(6,'(a)') 'Heun integration method selected'
#else
        write(6,'(a)') 'Heun integration method selected'
#endif
       case(2)
!-----------------------
! SIA
#ifdef CPP_MPI
        if (irank.eq.0) write(6,'(a)') 'SIA: Heun semi-implicit+projector integration method selected'
#else
        write(6,'(a)') 'SIA: Heun semi-implicit+projector integration method selected'
#endif
       Bini=0.0d0
       case(4)
!-----------------------
! SIB
#ifdef CPP_MPI
        if (irank.eq.0) write(6,'(a)') 'SIB: Heun semi-implicit+projector+NC integration method selected'
#else
        write(6,'(a)') 'SIB: Heun semi-implicit+projector+NC integration method selected'
#endif
       Bini=0.0d0
       case(1)
!-----------------------
! Euler
#ifdef CPP_MPI1
        if (irank.eq.0) write(6,'(a)') 'Euler integration method selected'
#else
        write(6,'(a)') 'Euler integration method selected'
        
#endif
       Bini=0.0d0
       case(5)
!-----------------------
! SIB with error correction
#ifdef CPP_MPI
        if (irank.eq.0) write(6,'(a)') 'SIB with error correction method selected'
#else
        write(6,'(a)') 'SIB with error correction method selected'
#endif
       Bini=0.0d0
       case(6)
!-----------------------
! SIB without temperature
#ifdef CPP_MPI
        if (irank.eq.0) write(6,'(a)') 'SIB: Heun semi-implicit+projector+NC integration method selected (T=0)'
#else
        write(6,'(a)') 'SIB: Heun semi-implicit+projector+NC integration method selected (T=0)'
#endif
       Bini=0.0d0
       case default
      end select

      iomp=1
      iomp=maxloc(htor)
      maxh=maxval(htor)
      if (maxh.lt.1.0d-8) maxh=1.0d0
      check=0.0d0
      check1=0.0d0
      check2=0.0d0
      check3=0.0d0
      Eold=100.0d0
      k=0
      h=0
      l=0
      kt=ktini
      step=0.0d0
      steptor=0.0d0
      stepadia=0.0d0
      stepsttor=0.0d0
      steptemp=0.0d0
      real_time=0.0d0
      Einitial=0.0d0
      h_int=h_ext
      said_it_once=.False.

      stmtorquebp=storque
      security=0.0d0
#ifdef CPP_MPI
      if (irank.eq.0) then
#endif
      B=Bexch(iomp(1),iomp(2),iomp(3),1,spin,shape_spin,indexNN,shape_index,masque,shape_masque,tableNN,shape_tableNN) &
    &  +BZ(iomp(1),iomp(2),iomp(3),1,spin,shape_spin,h_int) &
    &  +Bani(iomp(1),iomp(2),iomp(3),1,EA,spin,shape_spin,masque,shape_masque) &
    &  +Efield(iomp(1),iomp(2),iomp(3),1,spin,shape_spin,indexNN,shape_index,masque,shape_masque,tableNN,shape_tableNN)
      if (i_dm) B=B+Bdm(iomp(1),iomp(2),iomp(3),1,spin,shape_spin,indexNN,shape_index,masque,shape_masque,tableNN,shape_tableNN)
      if (i_biq) B=B+Bbiqd(iomp(1),iomp(2),iomp(3),1,spin,shape_spin,indexNN,shape_index,masque,shape_masque,tableNN,shape_tableNN)
      if (i_four) B=B+Bfour(iomp(1),iomp(2),iomp(3),1,spin,shape_spin,masque,shape_masque)

      step=cross(B,spin(4:6,iomp(1),iomp(2),iomp(3),1))
      if (i_torque) then
          write(6,'(a)') 'you are now using the SHE+STT'
          write(6,'(a)') 'field like STT proportional to torque*(1-damping*damping_torque)'
          write(6,'(a)') 'Anti-field like STT proportional to torque*(damping+dampin_torque)'
          dumy(1:3)=cross(spin(4:6,iomp(1),iomp(2),iomp(3),1),Ipol)
          steptor=torque_FL*(1.0d0-damping*torque_AFL)*dumy(1:3)+ &
    &  torque_FL*(damping+torque_AFL)*cross(spin(4:6,iomp(1),iomp(2),iomp(3),1),dumy(1:3))
      endif
      if ((abs(nonadia).lt.1.0d-8).and.(abs(adia).lt.1.0d-8)) then
          write(6,'(a)') 'you are now using the STT for in plane current'
          stepadia=cross(-spin(4:6,tableNN(1,1,iomp(1),iomp(2),iomp(3),1),tableNN(2,1,iomp(1),iomp(2),iomp(3),1), &
       tableNN(3,1,iomp(1),iomp(2),iomp(3),1),1),spin(4:6,iomp(1),iomp(2),iomp(3),1))
      endif
      if (stmtorque) stepsttor=cross(spin(4:6,iomp(1),iomp(2),iomp(3),1),Ipol*htor(iomp(1),iomp(2),iomp(3)))
      if (stmtemp) then
       steptemp=cross(spin(4:6,iomp(1),iomp(2),iomp(3),1),(/randist(kt),randist(kt),randist(kt)/)) &
     & *htor(iomp(1),iomp(2),iomp(3))/maxh
      else
       steptemp=cross(spin(4:6,iomp(1),iomp(2),iomp(3),1),(/randist(kt),randist(kt),randist(kt)/))
      endif

       ds=timestep*(step+damping*cross(step,spin(4:6,iomp(1),iomp(2),iomp(3),1))+torque_FL* &
     &  cross(steptor,spin(4:6,iomp(1),iomp(2),iomp(3),1))+adia*&
     &  cross(spin(4:6,iomp(1),iomp(2),iomp(3),1),stepadia)-nonadia*stepadia   &
     &  +storque*cross(stepsttor,spin(4:6,iomp(1),iomp(2),iomp(3),1)))/hbar

      write(6,'(a)') 'The LLG equations are used'
      write(6,'(a)') ' '
      write(6,'(a,3I5,3x,a,I9)') 'orders of magnitude for spin ', iomp, 'milieu', product(shape_spin(2:4))/2
      write(6,'(a,f12.8,a,f9.4,a)') 'Beff ',norm(B),', damping ', &
        damping*100.0,'% of Beff'
      if (torque_FL.ne.0.0d0) write(6,'(2(a,f12.8,2x))') 'torque ', torque_FL, 'norm ', norm(steptor)
      if (stmtorque) write(6,'(2(a,f12.8,2x))') 'SPSTM-torque ',storque, 'norm ', norm(stepsttor)
      write(6,'(a,2x,f12.8)') 'ds', norm(ds)
      if (kT.ne.0.0d0) write(6,'(a,2x,f12.6,2x,a,2x,f14.10,2x,a,f12.9)') 'Temperature', kT/k_B &
     &   ,'amplitude', dsqrt(2.0d0*damping*kT), 'example dsT', dsqrt(2.0d0*damping*kT)*norm(steptemp)
      if (rampe) write(6,'(a,2x,f12.6,2x,a,2x,f12.7,2x,a,I6,2x,a,I6)') 'step rampe', temps/k_B &
     &   , 'Tf=', kTfin/k_B, 'time step', times, 'time start', tstart
      if (hsweep) write(6,'(a,2x,3(f8.4,2x),a,2x,3(f8.4,2x),a,I6,2x,a,I6,2x,a,3(f8.4,2x))') 'Hini', (H_int(i),i=1,3) &
     &   , 'Hstep', (hstep(i),i=1,3), 'time step', htimes, 'time start', hstart, 'Hfin', (hfin(i),i=1,3)

#ifdef CPP_MPI
      endif
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! beginning of the
      do j=1,duration
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifndef CPP_BRUTDIP
      if (i_dip) then
#ifdef CPP_OPENMP
!$OMP parallel private(i_x,i_y,i_z) default(shared)
#endif
       do i_z=1,shape_spin(4)
        do i_y=1,shape_spin(3)
         do i_x=1,shape_spin(2)
        mmatrix(:,i_x,i_y,i_z)=spin(4:6,i_x,i_y,i_z,1)
         enddo
        enddo
       enddo
#ifdef CPP_OPENMP
!$OMP end parallel
#endif
      endif
#endif
! send the lattice to the processors
!------------------------------
!        call dyna_split(i_dip,MPI_WORKING_WORLD,irank)
!       call MPI_BCAST(Spin(4:7,:,:,:,:),product(dim_lat)*count(motif%i_m)*4,MPI_REAL8,0, &
!     &     MPI_WORKING_WORLD,ierr)
!       write(*,*) '-----------------------'
!       write(*,*) j
!       write(*,*) '-----------------------'

       spinafter=0.0d0
       spinafter(4,:,:,:,:)=spin(7,1,1,1,1)
       qeuler=0.0d0
       vx=0.0d0
       vy=0.0d0
       vz=0.0d0
       Mx=0.0d0
       My=0.0d0
       Mz=0.0d0
       Edy=0.0d0
       Mdy=0.0d0
       vortex=0.0d0
       test_torque=0.0d0
       ave_torque=0.0d0
       l=l+1
       h=h+1

       if (((j.lt.ti).or.(j.gt.tf)).and.(marche)) then
        storque=0.0d0
        elseif (((j.gt.ti).and.(j.lt.tf)).and.(marche)) then
        storque=stmtorquebp
       endif

       if ((l.gt.times).and.(rampe).and.(j.gt.tstart)) then
        kT=kT+temps
        l=1
        write(6,'(a,2x,f16.6)') 'Final Temp', check(1)/check(2)/2.0d0/k_B
        write(6,'(a,f10.5)') 'T=',kT/k_B
        check=0.0d0
       endif
       if ((h.gt.htimes).and.(hsweep).and.(norm((H_int-Hfin)).gt.1.0d-6).and. &
          (j.gt.hstart)) then
        H_int=H_int+hstep
        h=1
#ifdef CPP_MPI
        if (irank.eq.0) write(6,'(a,2x,3f8.4)') 'applied field', (H_int(i),i=1,3)
#else
        write(6,'(a,2x,3f8.4)') 'applied field', (H_int(i),i=1,3)
#endif
       endif

! different integration types
!-----------------------------------------------
! Euler integration scheme
!-----------------------------------------------
      select case (integtype)
       case (1)

#ifdef CPP_OPENMP
!$OMP parallel private(i_x,i_y,i_z,i_m,Beff) default(shared) reduction(+:check1,check2,check3)
#endif
       test_torque=0.0d0
       do i_m=1,Mstop
        do i_z=Zstart,Zstop
         do i_y=Ystart,Ystop
          do i_x=Xstart,Xstop

       call calculate_Beff(i_DM,i_four,i_biq,i_dip,EA,i_x,i_y,i_z,i_m,Beff, &
          &  spin,shape_spin,indexNN,shape_index,masque,shape_masque,tableNN,shape_tableNN,h_int)

      spinafter(1:3,i_x,i_y,i_z,i_m)=integrate(timestep,Beff,kt,stmtemp,maxh,i_x,i_y,i_z &
     & ,i_m,damping,Ipol,torque_FL,torque_AFL,adia,nonadia,storque,i_torque,stmtorque,spin,shape_spin,tableNN)

! the temperature is checked with 1 temperature step before
!!! check temperature
       call update_temp_measure(check1,check2,spinafter(1:3,i_x,i_y,i_z,i_m),Beff)
       if (norm_cross(spinafter(1:3,i_x,i_y,i_z,i_m),Beff).gt.test_torque) test_torque=norm_cross(spinafter(1:3,i_x,i_y,i_z,i_m),Beff)
!!! end check

          enddo
         enddo
        enddo
       enddo

#ifdef CPP_OPENMP
!$OMP end parallel
#endif
       check(1)=check(1)+check1
       check(2)=check(2)+check2
       real_time=real_time+timestep
#ifdef CPP_MPI
       trans(1)=test_torque
       call MPI_REDUCE(trans(1),test_torque,1,MPI_REAL8,MPI_SUM,0,MPI_COMM,ierr)
#endif
       if (j.eq.1) check3=test_torque
!-----------------------------------------------
! Heun integration scheme
!-----------------------------------------------
       case (3)

#ifdef CPP_OPENMP
!$OMP parallel private(i_x,i_y,i_z,i_m) default(shared)
#endif
       do i_m=1,Mstop
        do i_z=Zstart,Zstop
         do i_y=Ystart,Ystop
          do i_x=Xstart,Xstop

          spinini(:,i_x,i_y,i_z,i_m)=spin(4:7,i_x,i_y,i_z,i_m)
          enddo
         enddo
        enddo
       enddo
#ifdef CPP_OPENMP
!$OMP end parallel

!$OMP parallel private(i_x,i_y,i_z,i_m) default(shared)
#endif
       do i_m=1,Mstop
        do i_z=Zstart,Zstop
         do i_y=Ystart,Ystop
          do i_x=Xstart,Xstop

       call calculate_Beff(i_DM,i_four,i_biq,i_dip,EA,i_x,i_y,i_z,i_m,Bini(:,i_x,i_y,i_z,i_m), &
  &      spin,shape_spin,indexNN,shape_index,masque,shape_masque,tableNN,shape_tableNN,h_int)

       spinafter(1:3,i_x,i_y,i_z,i_m)=integrate(timestep,Bini(:,i_x,i_y,i_z,i_m),kt,stmtemp, &
     & maxh,i_x,i_y,i_z,i_m,damping,Ipol,torque_FL,torque_AFL,adia,nonadia,storque,i_torque,stmtorque,spin,shape_spin,tableNN)

       spinafter(4,i_x,i_y,i_z,i_m)=spinini(4,i_x,i_y,i_z,i_m)

          enddo
         enddo
        enddo
       enddo
#ifdef CPP_OPENMP
!$OMP end parallel
#endif

#ifdef CPP_MPI
! gather to take into account the possible change of size of spinafter
!------------------------------
       if (i_ghost) then
        call rebuild_mat(spinafter,N_site_comm*4,spin)
       else
        spin(4:7,:,:,:,:)=spinafter(1:4,:,:,:,:)
       endif

#else
! transfer of the predicator in the spin for the calculation of Beff
#ifdef CPP_OPENMP
!$OMP parallel private(i_x,i_y,i_z,i_m) default(shared)
#endif
       do i_m=1,Mstop
        do i_z=1,shape_spin(4)
         do i_y=1,shape_spin(3)
          do i_x=1,shape_spin(2)
        spin(4:7,i_x,i_y,i_z,i_m)=spinafter(:,i_x,i_y,i_z,i_m)
          enddo
         enddo
        enddo
       enddo
#ifdef CPP_OPENMP
!$OMP end parallel
#endif
#endif

#ifdef CPP_OPENMP
!$OMP parallel private(i_x,i_y,i_z,i_m,Beff) default(shared) reduction(+:check1,check2,check3)
#endif
       test_torque=0.0d0

       do i_m=1,Mstop
        do i_z=Zstart,Zstop
         do i_y=Ystart,Ystop
          do i_x=Xstart,Xstop

        call calculate_Beff(i_DM,i_four,i_biq,i_dip,EA,i_x,i_y,i_z,i_m,Beff, &
  &      spin,shape_spin,indexNN,shape_index,masque,shape_masque,tableNN,shape_tableNN,h_int)

        Beff=(Beff+Bini(:,i_x,i_y,i_z,i_m))/2.0d0

        spinafter(1:3,i_x,i_y,i_z,i_m)=integrate(timestep,Beff,kt,stmtemp,maxh,i_x,i_y,i_z &
      & ,i_m,damping,Ipol,torque_FL,torque_AFL,adia,nonadia,storque,i_torque,stmtorque,spinini,shape_spinini,tableNN)

        spinafter(4,i_x,i_y,i_z,i_m)=spinini(4,i_x,i_y,i_z,i_m)

! the temperature is checked with 1 temperature step before
!!! check temperature
       call update_temp_measure(check1,check2,spinafter(1:3,i_x,i_y,i_z,i_m),Beff)
       if (norm_cross(spinafter(1:3,i_x,i_y,i_z,i_m),Beff).gt.test_torque) test_torque=norm_cross(spinafter(1:3,i_x,i_y,i_z,i_m),Beff)
!!! end check

          enddo
         enddo
        enddo
       enddo
#ifdef CPP_OPENMP
!$OMP end parallel
#endif
       check(1)=check(1)+check1
       check(2)=check(2)+check2
       real_time=real_time+timestep
#ifdef CPP_MPI
       trans(1)=test_torque
       call MPI_REDUCE(trans(1),test_torque,1,MPI_REAL8,MPI_SUM,0,MPI_COMM,ierr)
#endif
       if (j.eq.1) check3=test_torque
!-----------------------------------------------
! SIA and IMP integration scheme
!-----------------------------------------------
       case (2)
! the position 1 of the predicator is the spins at time 0
#ifdef CPP_OPENMP
!$OMP parallel private(i_x,i_y,i_z,i_m) default(shared)
#endif

       do i_m=1,Mstop
        do i_z=1,shape_spin(4)
         do i_y=1,shape_spin(3)
          do i_x=1,shape_spin(2)
       spinini(:,i_x,i_y,i_z,i_m)=spin(4:7,i_x,i_y,i_z,i_m)
          enddo
         enddo
        enddo
       enddo
#ifdef CPP_OPENMP
!$OMP end parallel
#endif

#ifdef CPP_OPENMP
!$OMP parallel private(i_x,i_y,i_z,i_m,Beff) default(shared)
#endif
       do i_m=1,Mstop
        do i_z=Zstart,Zstop
         do i_y=Ystart,Ystop
          do i_x=Xstart,Xstop

       call calculate_Beff(i_DM,i_four,i_biq,i_dip,EA,i_x,i_y,i_z,i_m,Beff, &
  &      spin,shape_spin,indexNN,shape_index,masque,shape_masque,tableNN,shape_tableNN,h_int)

        spinafter(1:3,i_x,i_y,i_z,i_m)=(integrate(timestep,spin(4:6,i_x,i_y,i_z,i_m),Beff,kt,damping &
     & ,stmtemp,state,i_torque,stmtorque,torque_FL,torque_AFL,adia,nonadia,storque,maxh,Ipol,i_x,i_y,i_z,i_m,spin)+ &
     & spinini(1:3,i_x,i_y,i_z,i_m))/2.0d0

        spinafter(4,i_x,i_y,i_z,i_m)=spinini(4,i_x,i_y,i_z,i_m)

          enddo
         enddo
        enddo
       enddo

#ifdef CPP_OPENMP
!$OMP end parallel
#endif

#ifdef CPP_MPI
! reduce
!------------------------------
       if (i_ghost) call rebuild_mat(spinafter,N_site_comm*4,spin)

#else
! transfer of the predicator in the spin for the calculation of Beff
        spin(4:6,:,:,:,:)=spinafter(1:3,:,:,:,:)
        spin(7,:,:,:,:)=spinini(4,:,:,:,:)
#endif

       test_torque=0.0d0

#ifdef CPP_OPENMP
!$OMP parallel private(i_x,i_y,i_z,i_m,Beff) default(shared) reduction(max:test_torque)
#endif
       do i_m=1,Mstop
        do i_z=Zstart,Zstop
         do i_y=Ystart,Ystop
          do i_x=Xstart,Xstop

        call calculate_Beff(i_DM,i_four,i_biq,i_dip,EA,i_x,i_y,i_z,i_m,Beff, &
     &    spin,shape_spin,indexNN,shape_index,masque,shape_masque,tableNN,shape_tableNN,h_int)

        spinafter(1:3,i_x,i_y,i_z,i_m)=integrate(timestep,spinini(1:3,i_x,i_y,i_z,i_m),Beff,kt,damping &
     &   ,stmtemp,state,i_torque,stmtorque,torque_FL,torque_AFL,adia,nonadia,storque,maxh,check,Ipol,i_x,i_y,i_z,i_m,spin)

          spinafter(4,i_x,i_y,i_z,i_m)=spinini(4,i_x,i_y,i_z,i_m)

          if (norm_cross(spinafter(1:3,i_x,i_y,i_z,i_m),Beff).gt.test_torque) test_torque=norm_cross(spinafter(1:3,i_x,i_y,i_z,i_m),Beff)

          enddo
         enddo
        enddo
       enddo

#ifdef CPP_OPENMP
!$OMP end parallel
#endif

       real_time=real_time+timestep
#ifdef CPP_MPI
       trans(1)=test_torque
       call MPI_REDUCE(trans(1),test_torque,1,MPI_REAL8,MPI_MAX,0,MPI_COMM,ierr)
#endif
       if (j.eq.1) check3=test_torque
!-----------------------------------------------
! SIB and IMP integration scheme
!-----------------------------------------------
       case (4)
! the position 1 of the predicator is the spins at time 0

#ifdef CPP_OPENMP
!$OMP parallel private(i_x,i_y,i_z,i_m) default(shared)
#endif

       do i_m=1,Mstop
        do i_z=1,shape_spin(4)
         do i_y=1,shape_spin(3)
          do i_x=1,shape_spin(2)
       spinini(:,i_x,i_y,i_z,i_m)=spin(4:7,i_x,i_y,i_z,i_m)
          enddo
         enddo
        enddo
       enddo
#ifdef CPP_OPENMP
!$OMP end parallel
#endif

#ifdef CPP_OPENMP
!$OMP parallel private(i_x,i_y,i_z,i_m,Beff) default(shared)
#endif

       do i_m=1,Mstop
        do i_z=Zstart,Zstop
         do i_y=Ystart,Ystop
          do i_x=Xstart,Xstop

       call calculate_Beff(i_DM,i_four,i_biq,i_dip,EA,i_x,i_y,i_z,i_m,Beff, &
     &  spin,shape_spin,indexNN,shape_index,masque,shape_masque,tableNN,shape_tableNN,h_int)

        spinafter(1:3,i_x,i_y,i_z,i_m)=(integrate(timestep,spin(4:6,i_x,i_y,i_z,i_m),Beff,kt,damping &
     & ,stmtemp,state,i_torque,stmtorque,torque_FL,torque_AFL,adia,nonadia,storque,maxh,check,Ipol,i_x,i_y,i_z,i_m,spin)+ &
     &  spinini(1:3,i_x,i_y,i_z,i_m))/2.0d0

          spinafter(4,i_x,i_y,i_z,i_m)=spinini(4,i_x,i_y,i_z,i_m)
          enddo
         enddo
        enddo
       enddo

#ifdef CPP_OPENMP
!$OMP end parallel
#endif

#ifdef CPP_MPI
! reduce
!------------------------------
       if (i_ghost) call rebuild_mat(spinafter,N_site_comm*4,spin)
#else

! transfer of the predicator in the spin for the calculation of Beff
        spin(4:6,:,:,:,:)=spinafter(1:3,:,:,:,:)
        spin(7,:,:,:,:)=spinini(4,:,:,:,:)
#endif
       test_torque=0.0d0
#ifdef CPP_OPENMP
!$OMP parallel private(i_x,i_y,i_z,i_m,Beff) default(shared) reduction(max:test_torque)
#endif

       do i_m=1,Mstop
        do i_z=Zstart,Zstop
         do i_y=Ystart,Ystop
          do i_x=Xstart,Xstop

        call calculate_Beff(i_DM,i_four,i_biq,i_dip,EA,i_x,i_y,i_z,i_m,Beff, &
     &  spin,shape_spin,indexNN,shape_index,masque,shape_masque,tableNN,shape_tableNN,h_int)

        spinafter(1:3,i_x,i_y,i_z,i_m)=integrate(timestep,spinini(1:3,i_x,i_y,i_z,i_m),Beff,kt,damping &
     & ,stmtemp,state,i_torque,stmtorque,torque_FL,torque_AFL,adia,nonadia,storque,maxh,check,Ipol,i_x,i_y,i_z,i_m,spin)


          spinafter(4,i_x,i_y,i_z,i_m)=spinini(4,i_x,i_y,i_z,i_m)

          if (norm_cross(spinafter(1:3,i_x,i_y,i_z,i_m),Beff).gt.test_torque) test_torque=norm_cross(spinafter(1:3,i_x,i_y,i_z,i_m),Beff)

          enddo
         enddo
        enddo
       enddo
       real_time=real_time+timestep

#ifdef CPP_OPENMP
!$OMP end parallel
#endif

#ifdef CPP_MPI
       trans(1)=test_torque
       call MPI_REDUCE(trans(1),test_torque,1,MPI_REAL8,MPI_MAX,0,MPI_COMM,ierr)
#endif
       if (j.eq.1) check3=test_torque

!-----------------------------------------------
! SIB without temperature and with error control
!-----------------------------------------------
       case (5)

! save the initial spin state
       spinini=spin(4:7,:,:,:,:)

       call error_correction_SD(timestep,max_error,real_time,h_int, &
     & damping,i_torque,stmtorque,torque_FL,torque_AFL,adia,nonadia,storque,maxh,Ipol,N_site_comm,check3,j &
     & ,spin,shape_spin,indexNN,shape_index,masque,shape_masque,tableNN,shape_tableNN, &
     & i_biq,i_dm,i_four,i_dip,EA)

       if (timestep.lt.timestep_ini/20.0d0) then
        write(6,'(a)') 'the time step has decreased by a factor of 20'
        write(6,'(a)') 'the simulation might be very wrong'
        stop
       endif

!-----------------------------------------------
! SIB and IMP integration scheme without temperature
!-----------------------------------------------
       case (6)
! the position 1 of the predicator is the spins at time 0

#ifdef CPP_OPENMP
!$OMP parallel private(i_x,i_y,i_z,i_m) default(shared)
#endif

       do i_m=1,Mstop
        do i_z=1,shape_spin(4)
         do i_y=1,shape_spin(3)
          do i_x=1,shape_spin(2)
       spinini(:,i_x,i_y,i_z,i_m)=spin(4:7,i_x,i_y,i_z,i_m)
          enddo
         enddo
        enddo
       enddo
#ifdef CPP_OPENMP
!$OMP end parallel
#endif

#ifdef CPP_OPENMP
!$OMP parallel private(i_x,i_y,i_z,i_m,Beff) default(shared)
#endif
       do i_m=1,Mstop
        do i_z=Zstart,Zstop
         do i_y=Ystart,Ystop
          do i_x=Xstart,Xstop

       call calculate_Beff(i_DM,i_four,i_biq,i_dip,EA,i_x,i_y,i_z,i_m,Beff, &
     &   spin,shape_spin,indexNN,shape_index,masque,shape_masque,tableNN,shape_tableNN,h_int)

! check if there is any thing to do
        dum_norm=sqrt(Beff(1)**2+Beff(2)**2+Beff(3)**2)
        if ((dum_norm.lt.1.0d-8).or.((abs(dot_product(Beff,spin(4:6,i_x,i_y,i_z,i_m))/dum_norm-1.0d0)).lt.1.0d-10)) then
         spinafter(:,i_x,i_y,i_z,i_m)=spinini(:,i_x,i_y,i_z,i_m)
         cycle
        endif

        spinafter(1:3,i_x,i_y,i_z,i_m)=(integrate(timestep,spin(4:6,i_x,i_y,i_z,i_m),Beff,damping &
     & ,i_torque,stmtorque,torque_FL,torque_AFL,adia,nonadia,storque,maxh,Ipol,i_x,i_y,i_z,i_m,spin)+ &
     &  spinini(1:3,i_x,i_y,i_z,i_m))/2.0d0

!        if ((i_x.eq.1).and.(i_y.eq.1)) then
!         write(*,*) spinafter(1:3,i_x,i_y,i_z,i_m), norm(spinafter(1:3,i_x,i_y,i_z,i_m))
!         write(*,*) Exchange(i_x,i_y,i_z,i_m,spin,tableNN,masque,indexNN),Zeeman(i_x,i_y,i_z,i_m,spin,masque),anisotropy(i_x,i_y,i_z,i_m,EA,spin,masque)
!         write(*,*) Exchange(i_x,i_y,i_z,i_m,spin,tableNN,masque,indexNN)+Zeeman(i_x,i_y,i_z,i_m,spin,masque)+anisotropy(i_x,i_y,i_z,i_m,EA,spin,masque)
!
!         spin(4:6,i_x,i_y,i_z,i_m)=spinafter(1:3,i_x,i_y,i_z,i_m)
!         write(*,*) Exchange(i_x,i_y,i_z,i_m,spin,tableNN,masque,indexNN),Zeeman(i_x,i_y,i_z,i_m,spin,masque),anisotropy(i_x,i_y,i_z,i_m,EA,spin,masque)
!         write(*,*) Exchange(i_x,i_y,i_z,i_m,spin,tableNN,masque,indexNN)+Zeeman(i_x,i_y,i_z,i_m,spin,masque)+anisotropy(i_x,i_y,i_z,i_m,EA,spin,masque)
!
!         pause
!        endif

          spinafter(4,i_x,i_y,i_z,i_m)=spinini(4,i_x,i_y,i_z,i_m)
          enddo
         enddo
        enddo
       enddo

#ifdef CPP_OPENMP
!$OMP end parallel
#endif

#ifdef CPP_MPI
! reduce
!------------------------------
       if (i_ghost) call rebuild_mat(spinafter,N_site_comm*4,spin)
#else

! transfer of the predicator in the spin for the calculation of Beff
        spin(4:7,:,:,:,:)=spinafter(1:4,:,:,:,:)
#endif
       test_torque=0.0d0
#ifdef CPP_OPENMP
!$OMP parallel private(i_x,i_y,i_z,i_m,Beff) default(shared) reduction(max:test_torque) reduction(+:ave_torque)
#endif
       do i_m=1,Mstop
        do i_z=Zstart,Zstop
         do i_y=Ystart,Ystop
          do i_x=Xstart,Xstop

        call calculate_Beff(i_DM,i_four,i_biq,i_dip,EA,i_x,i_y,i_z,i_m,Beff, &
     &   spin,shape_spin,indexNN,shape_index,masque,shape_masque,tableNN,shape_tableNN,h_int)

! check if there is any thing to do
        dum_norm=sqrt(Beff(1)**2+Beff(2)**2+Beff(3)**2)
        if ((dum_norm.lt.1.0d-8).or.((abs(dot_product(Beff,spin(4:6,i_x,i_y,i_z,i_m))/dum_norm-1.0d0)).lt.1.0d-10)) then
         spinafter(:,i_x,i_y,i_z,i_m)=spinini(:,i_x,i_y,i_z,i_m)
         cycle
        endif


        spinafter(1:3,i_x,i_y,i_z,i_m)=integrate(timestep,spinini(1:3,i_x,i_y,i_z,i_m),Beff,damping &
     & ,i_torque,stmtorque,torque_FL,torque_AFL,adia,nonadia,storque,maxh,Ipol,i_x,i_y,i_z,i_m,spin)


          spinafter(4,i_x,i_y,i_z,i_m)=spinini(4,i_x,i_y,i_z,i_m)

          ave_torque=ave_torque+norm_cross(spinafter(1:3,i_x,i_y,i_z,i_m),Beff)
          if (norm_cross(spinafter(1:3,i_x,i_y,i_z,i_m),Beff).gt.test_torque) test_torque=norm_cross(spinafter(1:3,i_x,i_y,i_z,i_m),Beff)

          enddo
         enddo
        enddo
       enddo
       real_time=real_time+timestep
#ifdef CPP_OPENMP
!$OMP end parallel
#endif
#ifdef CPP_MPI
       trans(1)=test_torque
       call MPI_ALLREDUCE(trans(1),test_torque,1,MPI_REAL8,MPI_MAX,MPI_COMM,ierr)
       trans(1)=ave_torque
       call MPI_ALLREDUCE(trans(1),ave_torque,1,MPI_REAL8,MPI_SUM,MPI_COMM,ierr)
#endif
       if (j.eq.1) check3=test_torque

!-----------------------------------------------
       case default
       write(6,'(a)') 'no integration scheme chosen.'
       write(6,'(a)') 'check dyna.in file'
       stop
      end select

#ifdef CPP_MPI
! reduce
!------------------------------
       if (i_ghost)then
        call rebuild_mat(spinafter,N_site_comm*4,spin)
!        mpi_check=check
!        call MPI_ALLREDUCE(mpi_check,check,2,MPI_REAL8,MPI_SUM,MPI_COMM,ierr)
       else
        spin(4:7,:,:,:,:)=spinafter(1:4,:,:,:,:)
       endif

#else
#ifdef CPP_OPENMP
!$OMP parallel private(i_x,i_y,i_z,i_m) default(shared)
#endif
      do i_m=1,Mstop
       do i_z=1,shape_spin(4)
        do i_y=1,shape_spin(3)
         do i_x=1,shape_spin(2)
       spin(4:6,i_x,i_y,i_z,i_m)=spinafter(1:3,i_x,i_y,i_z,i_m)
         enddo
        enddo
       enddo
      enddo
#ifdef CPP_OPENMP
!$OMP end parallel
#endif
#endif

! calculate energy

#ifdef CPP_OPENMP
!$OMP parallel do private(i_x,i_y,i_z,i_m) default(shared) reduction(+:Edy,Mx,My,Mz,qeuler,vx,vy,vz)
#endif
      do i_m=1,Mstop
       do i_z=Zstart,Zstop
        do i_y=Ystart,Ystop
         do i_x=Xstart,Xstop

         if (masque(1,i_x,i_y,i_z).eq.0) cycle
         Edy=Edy+local_energy(Edy,i_DM,i_four,i_biq,i_dip,EA,i_x,i_y,i_z,i_m, &
              & spin,shape_spin,tableNN,shape_tableNN,masque,shape_masque,indexNN,shape_index,h_int)

         Mx=Mx+Spin(4,i_x,i_y,i_z,i_m)
         My=My+Spin(5,i_x,i_y,i_z,i_m)
         Mz=Mz+Spin(6,i_x,i_y,i_z,i_m)

         dumy=sd_charge(i_x,i_y,i_z,i_m,spin)

         qeuler=qeuler+dumy(1)/pi(4.0d0)
         vx=vx+dumy(2)
         vy=vy+dumy(3)
         vz=vz+dumy(4)
         enddo
        enddo
       enddo
      enddo
#ifdef CPP_OPENMP
!$OMP end parallel do
#endif
      Mdy=(/Mx,My,Mz/)
      vortex=(/vx,vy,vz/)/3.0d0/dsqrt(3.0d0)

#ifdef CPP_MPI
       trans(1)=Edy
       call MPI_ALLREDUCE(trans(1),Edy,1,MPI_REAL8,MPI_SUM,MPI_COMM,ierr)
       trans(1:3)=Mdy
       call MPI_REDUCE(trans(1:3),Mdy,3,MPI_REAL8,MPI_SUM,0,MPI_COMM,ierr)
       trans(1)=qeuler
       call MPI_REDUCE(trans(1),qeuler,1,MPI_REAL8,MPI_SUM,0,MPI_COMM,ierr)
       trans(1:3)=vortex
       call MPI_REDUCE(trans(1:3),vortex,3,MPI_REAL8,MPI_SUM,0,MPI_COMM,ierr)
       trans(1:2)=check
       call MPI_REDUCE(trans(1:2),check,2,MPI_REAL8,MPI_SUM,0,MPI_COMM,ierr)

      if (j.eq.1) Einitial=Edy/N_cell
#endif
#ifdef CPP_MPI
      if ((i_Efield).and.(mod(j-1,gra_freq).eq.0).and.(irank.eq.0)) call Efield_sd(j/gra_freq,spin,shape_spin,tableNN,shape_tableNN,masque,indexNN,h_int,irank,start,isize,MPI_COMM)
#else
      if ((i_Efield).and.(mod(j-1,gra_freq).eq.0)) call Efield_sd(j/gra_freq,spin,shape_spin,tableNN,masque,indexNN,h_int)
#endif

#ifdef CPP_MPI
       if (irank.eq.0) then
#endif
      Edy=Edy/N_cell
      Mdy=Mdy/N_cell

       if (dabs(check(2)).gt.1.0d-8) call get_temp(security,check,kt)

       if (mod(j-1,gra_freq).eq.0) Write(7,'(18(E20.12E3,2x),E20.12E3)') real_time,Edy, &
     &   norm(Mdy),Mdy,norm(vortex),vortex,qeuler, &
     &   storque, torque_FL, kT/k_B,(security(i),i=1,2),H_int

      if ((gra_log).and.(mod(j-1,gra_freq).eq.0)) then
         call CreateSpinFile(j/gra_freq,spin,shape_spin)
         call WriteSpinAndCorrFile(j/gra_freq,spin,shape_spin)
         write(6,'(a,I10)')'wrote Spin configuration and povray file number',j
      endif

      if ((gra_topo).and.(mod(j-1,gra_freq).eq.0)) then
       if (size(world).eq.1) then
        Call topocharge_sd(j/gra_freq,spin(1:6,:,1,1,1))
       elseif ((size(world).eq.2).and.(shape_spin(5).eq.1)) then
        Call topocharge_sd(j/gra_freq,spin(1:6,:,:,1,1))
       elseif ((size(world).eq.2).and.(shape_spin(5).ne.1)) then
        Call topocharge_sd(j/gra_freq,spin(1:6,:,:,1,:))
       else
        Call topocharge_sd(j/gra_freq,spin(1:6,:,:,:,:))
       endif
      endif

      if ((Ffield).and.(mod(j-1,gra_freq).eq.0)) call field_sd(j/gra_freq,spin,shape_spin,indexNN,shape_index,masque,shape_masque,tableNN,shape_tableNN,h_ext)

! security in case of energy increase in SD and check for convergence
      if (((damping*(Edy-Eold).gt.1.0d-10).or.(damping*(Edy-Einitial).gt.1.0d-10)).and.(kt.lt.1.0d-10).and.(.not.said_it_once)) then
#ifdef CPP_MPI
        write(6,'(a)') 'WARNING: the total energy or torque is increasing for non zero damping'
        write(6,'(a)') 'this is not allowed by theory'
        write(6,'(a)') 'please reduce the time step'
#else
       write(6,'(a)') 'WARNING: the total energy or torque is increasing for non zero damping'
       write(6,'(a)') 'this is not allowed by theory'
       write(6,'(a)') 'please reduce the time step'
#endif
       said_it_once=.True.
      endif

      if (mod(j-1,Efreq).eq.0) write(8,'(I10,3x,3(E20.12E3,3x))') j,Edy,test_torque,ave_torque
      check3=test_torque
      Eold=Edy

#ifdef CPP_MPI
      endif
      call mpi_barrier(MPI_COMM,ierr)
#endif

!!!! convergence criteria
!      if (test_torque.lt.1.0d-8) then
!       call cpu_time(computation_time)
!       write(*,*) 'computation time:',computation_time,'seconds'
!       write(*,*) 'simulation is converged'
!       exit
!      endif
        
!!!!!!!!!!!!!!! end of a timestep
       enddo

#ifdef CPP_MPI
      if (irank.eq.0) then
#endif
      close(7)
      close(8)

#ifdef CPP_MPI
      endif
#endif
      if ((dabs(check(2)).gt.1.0d-8).and.(kt/k_B.gt.1.0d-5)) then
       write(6,'(a,2x,f16.6)') 'Final Temp (K)', check(1)/check(2)/2.0d0/k_B
       write(6,'(a,2x,f14.7)') 'Kinetic energy (meV)', (check(1)/check(2)/2.0d0-kT)/k_B*1000.0d0
      else
       write(6,'(a)') 'the temperature measurement is not possible'
      endif

      end subroutine spindynamics
