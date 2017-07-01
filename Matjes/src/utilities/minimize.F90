
      subroutine minimize(i_biq,i_dm,i_four,i_dip,gra_log,gra_freq,EA, &
          & spin,shape_spin,tableNN,shape_tableNN,masque,shape_masque,indexNN,shape_index,N_cell,h_ext)
      use m_eval_Beff
      use m_write_spin
      use m_createspinfile
      use m_solver, only : minimization
      use m_local_energy
      use m_vector, only : norm_cross,norm,calculate_damping
#ifdef CPP_OPENMP
      use omp_lib
#endif
      implicit none
      logical, intent(in) :: i_biq,i_dm,i_four,i_dip,gra_log
      real(kind=8), intent(in) :: EA(3)
      integer, intent(in) :: shape_index(2),shape_spin(5),shape_tableNN(6),shape_masque(4),N_cell,gra_freq
      real(kind=8), intent(inout) :: spin(shape_spin(1),shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5))
      integer, intent(in) :: tableNN(shape_tableNN(1),shape_tableNN(2),shape_tableNN(3),shape_tableNN(4),shape_tableNN(5),shape_tableNN(6))
      integer, intent(in) :: masque(shape_masque(1),shape_masque(2),shape_masque(3),shape_masque(4)),indexNN(shape_index(1),shape_index(2))
      real(kind=8), intent(in) :: h_ext(3)
! dummy variable
      real(kind=8) :: velocity(3,shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5))
      real(kind=8) :: predicator(3,shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5))
      real(kind=8) :: force(3,shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5))
      real(kind=8) :: dt,masse,F_eff(3),V_eff(3),dumy,force_norm,Energy,max_torque,test_torque,conv_torque,vmax,vtest,F_temp(3)
! the computation time
      integer :: N_minimization,i_min,fin
      integer :: i_x,i_y,i_z,i_m
      logical :: i_exist
      character(len=100) :: str
      character(len=50) :: fname
      character(len=10) :: dummy
      integer, parameter :: io=70
#ifdef CPP_OPENMP
      integer :: nthreads,ithread
#endif


      masse=1.0d0
      dt=0.1d0
      N_minimization=1000
      max_torque=0.0d0
      test_torque=0.0d0
      conv_torque=1.0d-6

! read the input file
      i_exist=.False.
      inquire (file='minimization.in',exist=i_exist)

      if (.not.i_exist) then
       inquire (file='GNEB.in',exist=i_exist)
          if (i_exist) then
              write(6,'(a)') 'the minimization parameters are taken from file GNEB.in'
              write(6,'(a)') 'default parameters are'
              write(6,'(a,f8.6,2x,a,f8.6,2x,a,f10.7,2x,a,I6)') 'dt=',dt,'masse=',masse,'max torque=',conv_torque,'max steps=',N_minimization
              fname='GNEB.in'
          else
              write(6,'(a)') 'The file minimization.in was deleted or you are in the wrong routine'
              STOP
          endif
      else
        fname='minimization.in'
      endif

      open(io,file=fname,form='formatted',status='old',action='read')
      rewind(io)
      do
       read (io,'(a)',iostat=fin) str
       if (fin /= 0) exit
       str= trim(adjustl(str))
       if (len_trim(str)==0) cycle
       if (str(1:1) == '#' ) cycle

       if ( str(1:5) == 'steps') then
        backspace(io)
        read(io,*) dummy,N_minimization
       endif
       if ( str(1:5) == 'masse') then
        backspace(io)
        read(io,*) dummy,masse
       endif
       if ( str(1:8) == 'timestep') then
        backspace(io)
        read(io,*) dummy,dt
       endif
       if ( str(1:20) == 'convergence_criteria') then
        backspace(io)
        read(io,*) dummy,conv_torque
       endif

      enddo
      close(io)

!!! check the input
      if (masse.eq.0.0d0) then
       write(6,'(a)') 'The mass should be different from 0'
       stop
      endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      velocity=0.0d0
      force=0.0d0
      predicator=0.0d0
      F_eff=0.0d0
      V_eff=0.0d0
      force_norm=0.0d0
      vmax=0.0d0
      vtest=0.0d0

#ifdef CPP_OPENMP
nthreads=omp_get_num_procs()
call omp_set_num_threads(nthreads)

write(6,'(a)') 'OMP parallelization selected'
write(6,'(a,2x,I3,2x,a)') 'I will calculate on',nthreads,'threads'

#endif

      do i_m=1,shape_spin(5)
       do i_z=1,shape_spin(4)
        do i_y=1,shape_spin(3)
         do i_x=1,shape_spin(2)

           call calculate_Beff(i_DM,i_four,i_biq,i_dip,EA,i_x,i_y,i_z,i_m,F_eff, &
             & spin,shape_spin,indexNN,shape_index,masque,shape_masque,tableNN,shape_tableNN,h_ext)

           force(:,i_x,i_y,i_z,i_m)=calculate_damping(spin(4:6,i_x,i_y,i_z,i_m),F_eff)

           call minimization(spin(4:6,i_x,i_y,i_z,i_m),force(:,i_x,i_y,i_z,i_m),predicator(:,i_x,i_y,i_z,i_m),dt**2,masse*2.0d0,3)

           test_torque=norm_cross(predicator(:,i_x,i_y,i_z,i_m),force(:,i_x,i_y,i_z,i_m))

           if (test_torque.gt.max_torque) max_torque=test_torque

         enddo
        enddo
       enddo
      enddo

      do i_m=1,shape_spin(5)
       do i_z=1,shape_spin(4)
        do i_y=1,shape_spin(3)
         do i_x=1,shape_spin(2)
          spin(4:6,i_x,i_y,i_z,i_m)=predicator(:,i_x,i_y,i_z,i_m)
         enddo
        enddo
       enddo
      enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! end of initialization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(*,*) N_minimization
      do i_min=1,N_minimization

       max_torque=0.0d0
       dumy=0.0d0
       force_norm=0.0d0
       vmax=0.0d0

#ifdef CPP_OPENMP
!$OMP parallel default(shared) private(ithread)
ithread=omp_get_thread_num()

!$OMP do private(F_eff,V_eff) reduction(+:dumy,force_norm) reduction(max:vmax) collapse(4)

#endif

       do i_m=1,shape_spin(5)
        do i_z=1,shape_spin(4)
         do i_y=1,shape_spin(3)
          do i_x=1,shape_spin(2)

!          write(*,*) i_x,i_y,ithread

           call calculate_Beff(i_DM,i_four,i_biq,i_dip,EA,i_x,i_y,i_z,i_m,F_eff, &
               &   spin,shape_spin,indexNN,shape_index,masque,shape_masque,tableNN,shape_tableNN,h_ext)

           F_temp=calculate_damping(spin(4:6,i_x,i_y,i_z,i_m),F_eff)

           call minimization(velocity(:,i_x,i_y,i_z,i_m),(force(:,i_x,i_y,i_z,i_m)+F_temp)/2.0d0,V_eff,dt,masse,3)

           F_eff=F_temp
           force(:,i_x,i_y,i_z,i_m)=F_eff
           velocity(:,i_x,i_y,i_z,i_m)=V_eff

           dumy=dumy+V_eff(1)*F_eff(1)+V_eff(2)*F_eff(2)+V_eff(3)*F_eff(3)
           force_norm=force_norm+F_eff(1)**2+F_eff(2)**2+F_eff(3)**2

          enddo
         enddo
        enddo
       enddo

#ifdef CPP_OPENMP
!$OMP end do
#endif

       if (dumy.gt.0.0d0) then
#ifdef CPP_OPENMP
!$OMP do collapse(4)
#endif
        do i_m=1,shape_spin(5)
         do i_z=1,shape_spin(4)
          do i_y=1,shape_spin(3)
           do i_x=1,shape_spin(2)

             velocity(:,i_x,i_y,i_z,i_m)=dumy*force(:,i_x,i_y,i_z,i_m)/force_norm

           enddo
          enddo
         enddo
        enddo
#ifdef CPP_OPENMP
!$OMP end do
#endif
       else
            velocity=0.0d0
       endif

#ifdef CPP_OPENMP
!$OMP do private(i_x,i_y,i_z,i_m,test_torque) reduction(+:max_torque) collapse(4)
#endif

       do i_m=1,shape_spin(5)
        do i_z=1,shape_spin(4)
         do i_y=1,shape_spin(3)
          do i_x=1,shape_spin(2)

           call minimization(spin(4:6,i_x,i_y,i_z,i_m),velocity(:,i_x,i_y,i_z,i_m),force(:,i_x,i_y,i_z,i_m),predicator(:,i_x,i_y,i_z,i_m),dt,masse,3)

           test_torque=sqrt(force(1,i_x,i_y,i_z,i_m)**2+force(2,i_x,i_y,i_z,i_m)**2+force(3,i_x,i_y,i_z,i_m)**2)

           vtest=velocity(1,i_x,i_y,i_z,i_m)**2+velocity(2,i_x,i_y,i_z,i_m)**2+velocity(3,i_x,i_y,i_z,i_m)**2

           if (vtest.gt.vmax) vmax=vtest

           if (test_torque.gt.max_torque) max_torque=test_torque

           spin(4:6,i_x,i_y,i_z,i_m)=predicator(:,i_x,i_y,i_z,i_m)

          enddo
         enddo
        enddo
       enddo

#ifdef CPP_OPENMP
!$OMP end do
#endif
       Energy=0.0d0

#ifdef CPP_OPENMP
!$OMP do private(i_x,i_y,i_z,i_m) reduction(+:Energy) collapse(4)
#endif

       do i_m=1,shape_spin(5)
        do i_z=1,shape_spin(4)
         do i_y=1,shape_spin(3)
          do i_x=1,shape_spin(2)

           Energy=Energy+local_energy(Energy,i_DM,i_four,i_biq,i_dip,EA,i_x,i_y,i_z,i_m, &
              &   spin,shape_spin,tableNN,shape_tableNN,masque,shape_masque,indexNN,shape_index,h_ext)

          enddo
         enddo
        enddo
       enddo

#ifdef CPP_OPENMP
!$OMP end do
!$OMP end parallel
#endif

       

       if ((gra_log).and.(mod(i_min-1,gra_freq).eq.0)) then
         write(6,'(/,a,2x,I10)') 'iteration',i_min
         write(6,'(a,2x,f14.11)') 'Energy of the system (eV/unit cell)',Energy/dble(N_cell)
         write(6,'(2(a,2x,f14.11,2x))') 'convergence criteria:',conv_torque,',Measured Torque:',max_torque
         write(6,'(a,2x,f14.11,/)') 'speed of displacements:',vmax
          call WriteSpinAndCorrFile(i_min/gra_freq,spin,shape_spin)
          call CreateSpinFile(i_min/gra_freq,spin,shape_spin)
       endif

       if (conv_torque.gt.max_torque) then
        write(6,'(a)') 'minimization converged'
        exit
       endif

      enddo ! number of minimization steps

      end subroutine
