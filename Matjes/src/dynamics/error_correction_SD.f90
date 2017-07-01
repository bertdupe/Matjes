      module m_error_correction_SD
      interface error_correction_SD
       module procedure error_correction_SD_SIB
      end interface error_correction_SD
      contains

      subroutine error_correction_SD_SIB(timestep,max_error,real_time,h_ext, &
     & damping,i_torque,stmtorque,torque_FL,torque_AFL,adia,nonadia,storque,maxh,Ipol,N_site_comm,check3,j, &
     & spin,shape_spin,indexNN,shape_index,masque,shape_masque,tableNN,shape_tableNN, &
     & i_biq,i_dm,i_four,i_dip,EA)
      use m_solver
      use m_vector, only : cross,norm,norm_cross
      use m_constants, only : pi,k_b,hbar
      use m_eval_Beff
#ifndef CPP_BRUTDIP
      use m_setup_dipole, only : mmatrix
#endif
#ifdef CPP_MPI
     use m_parameters, only : i_ghost
     use m_make_box, only : Xstart,Xstop,Ystart,Ystop,Zstart,Zstop
     use m_mpi_prop, only : MPI_COMM,irank,isize
     use m_reconstruct_mat
#else
      use m_parameters, only : kt,ktini,ktfin,gra_topo,gra_log,gra_freq
#endif
      implicit none
! inout part of the variables
      logical, intent(in) :: i_biq,i_dm,i_four,i_dip
      real(kind=8), intent(in) :: max_error,h_ext(3),EA(3)
      real(kind=8), intent(in) :: damping,torque_FL,torque_AFL,adia,nonadia,storque,maxh,Ipol(:)
      logical, intent(in) :: i_torque,stmtorque
      integer, intent(in) :: N_site_comm,j
      integer, intent(in) :: shape_index(2),shape_spin(5),shape_tableNN(6),shape_masque(4)
      integer, intent(in) :: tableNN(shape_tableNN(1),shape_tableNN(2),shape_tableNN(3),shape_tableNN(4),shape_tableNN(5),shape_tableNN(6))
      integer, intent(in) :: masque(shape_masque(1),shape_masque(2),shape_masque(3),shape_masque(4))
      integer, intent(in) :: indexNN(shape_index(1),shape_index(2))
      real(kind=8), intent(inout) :: spin(shape_spin(1),shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5))
      real(kind=8), intent(inout) :: timestep,real_time,check3
! internal variables
      real(kind=8) :: spinend(4,shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5))
      real(kind=8) :: spinini(4,shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5))
      real(kind=8) :: spinafter(4,shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5))
      real(kind=8) :: spintest(4,shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5))
      real(kind=8) :: tfin,dt,Beff(3),error,max_timestep,test,test_torque
      integer :: i_x,i_y,i_z,i_m,ierr
#ifndef CPP_MPI
       integer, parameter :: Xstart=1
       integer, parameter :: Ystart=1
       integer, parameter :: Zstart=1
       integer :: Xstop,Ystop,Zstop

#else
       include 'mpif.h'
#endif

#ifndef CPP_MPI
      Xstop=shape_spin(2)
      Ystop=shape_spin(3)
      Zstop=shape_spin(4)
#endif
! error correction loop

       max_timestep=timestep*2.0d0
       dt=timestep
       error=100.0d0
       do while ((error.gt.max_error).and.(dt.ge.max_timestep/4.0d0))




!!!!!!!!!!!!!!!!!!!!!!!!!
! make the first big step
!!!!!!!!!!!!!!!!!!!!!!!!!



       spinini=spin(4:7,:,:,:,:)

! the position 1 of the predicator is the spins at time 0
       test_torque=0.0d0
       do i_m=1,shape_spin(5)

       do i_z=Zstart,Zstop
        do i_y=Ystart,Ystop
         do i_x=Xstart,Xstop

       call calculate_Beff(i_DM,i_four,i_biq,i_dip,EA,i_x,i_y,i_z,i_m,Beff,spin,shape_spin,indexNN,shape_index,masque,shape_masque,tableNN,shape_tableNN,h_ext)

        spinafter(1:3,i_x,i_y,i_z,i_m)=(integrate(2.0d0*dt,spin(4:6,i_x,i_y,i_z,i_m),Beff,damping &
     &  ,i_torque,stmtorque,torque_FL,torque_AFL,adia,nonadia,storque,maxh,Ipol,i_x,i_y,i_z,i_m,spin)+ &
     &  spinini(1:3,i_x,i_y,i_z,i_m))/2.0d0

          spinafter(4,i_x,i_y,i_z,i_m)=spinini(4,i_x,i_y,i_z,i_m)
          enddo
         enddo
        enddo
       enddo

#ifdef CPP_MPI
! reduce
!------------------------------
       if (i_ghost) call rebuild_mat(spinafter,N_site_comm*4,spin)
#else

! transfer of the predicator in the spin for the calculation of Beff
        spinini(1:3,:,:,:,:)=spinafter(1:3,:,:,:,:)
#endif

       do i_m=1,shape_spin(5)
       do i_z=Zstart,Zstop
        do i_y=Ystart,Ystop
         do i_x=Xstart,Xstop

        call calculate_Beff(i_DM,i_four,i_biq,i_dip,EA,i_x,i_y,i_z,i_m,Beff,spinini,shape_spin,indexNN,shape_index,masque,shape_masque,tableNN,shape_tableNN,h_ext)

        spinafter(1:3,i_x,i_y,i_z,i_m)=integrate(2.0d0*dt,spinini(1:3,i_x,i_y,i_z,i_m),Beff, &
     &   damping,i_torque,stmtorque,torque_FL,torque_AFL,adia,nonadia,storque,maxh,Ipol,i_x,i_y,i_z,i_m,spinini)

          spinafter(4,i_x,i_y,i_z,i_m)=spinini(4,i_x,i_y,i_z,i_m)

          test_torque=test_torque+norm_cross(spinafter(1:3,i_x,i_y,i_z,i_m),Beff)

          enddo
         enddo
        enddo
       enddo

#ifdef CPP_MPI
! reduce
!------------------------------
       if (i_ghost) then
        call rebuild_mat(spinafter,N_site_comm*4,spinend)
        call mpi_allreduce(test_torque,test,1,MPI_REAL8,MPI_SUM,MPI_COMM,ierr)
        test_torque=test
       endif
#else
! transfer of the predicator in the spin for the calculation of Beff
       spinend(:,:,:,:,:)=spinafter(1:4,:,:,:,:)
#endif
!! check that the torque should decrease
       if (j.eq.1) check3=test_torque
       write(*,*) 'beginning'
       write(*,*) abs(test_torque-check3),dt
       if (abs(test_torque-check3).gt.max_error) then
        dt=dt/2.0d0
        timestep=dt/2.0d0
        cycle
       endif




!!!!!!!!!!!!!!!!!!!!!!!!!
! make two small steps afterwards
!!!!!!!!!!!!!!!!!!!!!!!!!





       tfin=2.0d0*dt
       test_torque=0.0d0
       do while (dt.le.tfin)

! the position 1 of the predicator is the spins at time 0
       do i_m=1,shape_spin(5)
       do i_z=Zstart,Zstop
        do i_y=Ystart,Ystop
         do i_x=Xstart,Xstop

       call calculate_Beff(i_DM,i_four,i_biq,i_dip,EA,i_x,i_y,i_z,i_m,Beff,spinini,shape_spin,indexNN,shape_index,masque,shape_masque,tableNN,shape_tableNN,h_ext)

        spinafter(1:3,i_x,i_y,i_z,i_m)=(integrate(dt,spinini(1:3,i_x,i_y,i_z,i_m),Beff,damping &
     & ,i_torque,stmtorque,torque_FL,torque_AFL,adia,nonadia,storque,maxh,Ipol,i_x,i_y,i_z,i_m,spin)+ &
     &  spinini(1:3,i_x,i_y,i_z,i_m))/2.0d0

          spinafter(4,i_x,i_y,i_z,i_m)=spinini(4,i_x,i_y,i_z,i_m)
          enddo
         enddo
        enddo
       enddo

#ifdef CPP_MPI
! reduce
!------------------------------
       if (i_ghost) call rebuild_mat(spinafter,N_site_comm*4,spin)
#else

! transfer of the predicator in the spin for the calculation of Beff
        spin(4:6,:,:,:,:)=spinafter(1:3,:,:,:,:)
        spin(7,:,:,:,:)=spinini(4,:,:,:,:)
#endif

       do i_m=1,shape_spin(5)
       do i_z=Zstart,Zstop
        do i_y=Ystart,Ystop
         do i_x=Xstart,Xstop

        call calculate_Beff(i_DM,i_four,i_biq,i_dip,EA,i_x,i_y,i_z,i_m,Beff,spin,shape_spin,indexNN,shape_index,masque,shape_masque,tableNN,shape_tableNN,h_ext)

        spinafter(1:3,i_x,i_y,i_z,i_m)=integrate(dt,spinini(1:3,i_x,i_y,i_z,i_m),Beff, &
     &   damping,i_torque,stmtorque,torque_FL,torque_AFL,adia,nonadia,storque,maxh,Ipol,i_x,i_y,i_z,i_m,spin)

          spinafter(4,i_x,i_y,i_z,i_m)=spinini(4,i_x,i_y,i_z,i_m)

        test_torque=test_torque+norm_cross(spinafter(1:3,i_x,i_y,i_z,i_m),Beff)
          enddo
         enddo
        enddo
       enddo

#ifdef CPP_MPI
! reduce
!------------------------------
       if (i_ghost) call rebuild_mat(spinafter,N_site_comm*4,spinini)
#else
! transfer of the predicator in the spin for the calculation of Beff
       spinini=spinafter
#endif
       dt=2.0d0*dt
       enddo ! enddo loop for small dt

#ifdef CPP_MPI
       call mpi_allreduce(test_torque,test,1,MPI_REAL8,MPI_SUM,MPI_COMM,ierr)
       test_torque=test
#endif

!!!!!!!!!!!!!!!!
! check that the torque decreases before anything
!!!!!!!!!!!!!!!!


       write(*,*) 'torque'
       write(*,*) abs(test_torque-check3),max_error

       if (abs(test_torque-check3).gt.max_error) then
        write(6,'(a)') 'too big timestep: the torque is always increasing'
        write(6,*) 'initial timestep',max_timestep/2.0d0,'smaller timestep',timestep
        dt=dt/2.0d0
        timestep=dt
        cycle
       endif




!!!!!!!!!!!!!!!!!!
! test if the 2 end points are the same
!!!!!!!!!!!!!!!!!!


       spintest=spinini

       error=norm(spintest(:,Xstart,Ystart,Zstart,1)-spinend(:,Xstart,Ystart,Zstart,1))

       do i_m=1,shape_spin(5)
       do i_z=Zstart,Zstop
        do i_y=Ystart,Ystop
         do i_x=Xstart,Xstop

          test=norm(spintest(:,i_x,i_y,i_z,i_m)-spinend(:,i_x,i_y,i_z,i_m))
          if (test.gt.error) error=test

          enddo
         enddo
        enddo
       enddo

#ifdef CPP_MPI
       call mpi_allreduce(error,test,1,MPI_REAL8,MPI_MAX,MPI_COMM,ierr)
       error=test
#endif

       write(*,*) 'error'
        write(*,*) error,max_error
       if (error.lt.max_error) then
        spinafter=spinend
        timestep=dt
        real_time=real_time+2.0d0*dt
        exit
       else
        dt=timestep/2.0d0
       endif

       if (dt.lt.max_timestep/4.0d0) then
        spinafter=spintest
        timestep=max_timestep/4.0d0
        real_time=real_time+2.0d0*timestep
        write(6,'(a)') 'WARNING: your timestep is too large for the error you require'
        exit
       endif

      enddo ! enddo of the big loop over the timestep

!!!!!!!!!!!!!!
! case where it did not work at all
!!!!!!!!!!!!!!

      if (abs(test_torque-check3).gt.max_error) then
        write(6,'(a)') 'too big timestep: the torque is always increasing'
        write(6,*) 'initial timestep',max_timestep/2.0d0,'smaller timestep',timestep*2.0d0
        stop
      endif

! update the torque
      check3=test_torque
      write(6,*) 'The time step is now',timestep
      pause
      end subroutine error_correction_SD_SIB

      end module m_error_correction_SD
