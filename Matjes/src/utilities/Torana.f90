       module m_Torana
       use m_eval_Beff
       use m_vector
       use m_constants, only : pi
       use m_solver
       implicit none

       contains


       subroutine Torana(spin,shape_spin,indexNN,shape_index,masque,shape_masque,tableNN,shape_tableNN, &
                  torque_AFL,torque_FL,damping,Ipol,EA,h_int, &
                  i_DM,i_four,i_biq,i_dip,i_torque,stmtorque, &
                  timestep,adia,nonadia,storque,maxh)
       integer, intent(in) :: shape_index(2),shape_spin(5),shape_tableNN(6),shape_masque(4)
       real(kind=8), intent(in) :: spin(shape_spin(1),shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5))
       integer, intent(in) :: tableNN(shape_tableNN(1),shape_tableNN(2),shape_tableNN(3),shape_tableNN(4),shape_tableNN(5),shape_tableNN(6))
       integer, intent(in) :: masque(shape_masque(1),shape_masque(2),shape_masque(3),shape_masque(4)),indexNN(shape_index(1),shape_index(2))
       real(kind=8), intent(in) :: torque_AFL,torque_FL,Ipol(3),EA(3),h_int(3),damping,timestep,adia,nonadia,maxh,storque
       logical, intent(in) :: i_DM,i_four,i_biq,i_dip,i_torque,stmtorque
! internal variable
       integer :: i_m,i_z,i_y,i_x
       real(kind=8) :: Bini(3,shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5))
       real(kind=8) :: spinafter(4,shape_tableNN(3),shape_tableNN(4),shape_tableNN(5),shape_tableNN(6))
       real(kind=8) :: spinprint(6,shape_tableNN(3),shape_tableNN(4),shape_tableNN(5),shape_tableNN(6))
! angle for the maximum deviation between B_eff and M
       real(kind=8) :: max_angle,min_angle
       real(kind=8) :: dummy,test_ang
! for the mesh
       real(kind=8) :: T_FL_i,T_DL_i,pas_DL,pas_FL,Beff(3),torque_int_AFL,torque_int_FL
       character(len=30) :: fname,fname2,toto
       integer :: n_print_spin,k,i


       write(6,'(a,2x,3f12.6)') 'Current polarization direction',Ipol
       write(6,'(a,2x,f12.6)') 'Value of the field-like Torque (in eV)',torque_FL
       write(6,'(a,2x,f12.6,/)') 'Value of the damping-like Torque (in eV)',torque_FL*torque_AFL

! initialization of the variables

      Bini=0.0d0
      max_angle=0.0d0
      min_angle=10.0d0
      dummy=0.0d0
      test_ang=0.0d0
      spinprint=0.0d0
      n_print_spin=0
      torque_int_AFL=0.0d0
      torque_int_FL=0.0d0
      Beff=0.0d0

! first one hase to calculate the effective magnetic fields created by all the neighbours

       do i_m=1,shape_spin(5)
        do i_z=1,shape_spin(4)
         do i_y=1,shape_spin(3)
          do i_x=1,shape_spin(2)

       call calculate_Beff(i_DM,i_four,i_biq,i_dip,EA,i_x,i_y,i_z,i_m,Bini(:,i_x,i_y,i_z,i_m), &
  &      spin,shape_spin,indexNN,shape_index,masque,shape_masque,tableNN,shape_tableNN,h_int)

          dummy=(spin(4,i_x,i_y,i_z,i_m)*Bini(1,i_x,i_y,i_z,i_m)+spin(5,i_x,i_y,i_z,i_m)*Bini(2,i_x,i_y,i_z,i_m)+spin(6,i_x,i_y,i_z,i_m)*Bini(3,i_x,i_y,i_z,i_m))/norm(Bini(:,i_x,i_y,i_z,i_m))
          test_ang=acos(dummy)*180/pi(1.0d0)

          if (test_ang.gt.max_angle) max_angle=test_ang
          if (test_ang.lt.min_angle) min_angle=test_ang

          enddo
         enddo
        enddo
       enddo

       write(6,'(a,4x,f12.6)') 'MAXImum deviation angle between B_eff and M in degree',max_angle
       write(6,'(a,4x,f12.6,/)') 'MINImum deviation angle between B_eff and M in degree',min_angle

! !!! If the equations of motions are changed then this must be changed too !!!

! value of the effective magnetic field that contains the the Polarization of the current - field-like
! torque_FL*(1.0d0-damping*torque_AFL)*Ipol-B_eff
       write(6,'(a,4x,f12.6,/)') 'Amplitude of the deviations of field-like Torque (in eV)', torque_FL*(1.0d0-damping*torque_AFL)
       write(6,'(a,4x,3f12.6,/)') 'Direction of the current Polarization', Ipol

! Now the damping like torque
! torque_FL*(torque_AFL+damping)*Ipol-damping*B_eff
       write(6,'(a,4x,f12.6,/)') 'Amplitude of the deviations of damping-like Torque (in eV)', torque_FL*(torque_AFL+damping)

! Let's have a mesh purely done on the torque field like

! first test the numbers of spins which must be printed

       do i_m=1,shape_spin(5)
        do i_z=1,shape_spin(4)
         do i_y=1,shape_spin(3)
          do i_x=1,shape_spin(2)

           dummy=(spin(4,i_x,i_y,i_z,i_m)*Bini(1,i_x,i_y,i_z,i_m)+spin(5,i_x,i_y,i_z,i_m)*Bini(2,i_x,i_y,i_z,i_m)+spin(6,i_x,i_y,i_z,i_m)*Bini(3,i_x,i_y,i_z,i_m))/norm(Bini(:,i_x,i_y,i_z,i_m))
           test_ang=acos(dummy)*180/pi(1.0d0)

           if (test_ang.gt.min_angle+1.0d-4) n_print_spin=n_print_spin+1

          enddo
         enddo
        enddo
       enddo

       write(6,'(a,4x,I10,/)') 'Number of spins which will be printed', n_print_spin

       T_FL_i=torque_FL
       T_DL_i=torque_AFL
       pas_DL=torque_AFL/10.0d0
       pas_FL=torque_FL/10.0d0

       do k=0,10
        torque_int_FL=pas_FL*real(k)

        do i_m=1,shape_spin(5)
         do i_z=1,shape_spin(4)
          do i_y=1,shape_spin(3)
           do i_x=1,shape_spin(2)

           Beff=Bini(:,i_x,i_y,i_z,i_m)

!           spinafter(1:3,i_x,i_y,i_z,i_m)=integrate(timestep,Beff,0.0,.False.,maxh,i_x,i_y,i_z &
!     & ,i_m,damping,Ipol,torque_int_FL,torque_AFL,adia,nonadia,storque,i_torque,stmtorque,spin,shape_spin,tableNN)

           dummy=(spin(4,i_x,i_y,i_z,i_m)*Bini(1,i_x,i_y,i_z,i_m)+spin(5,i_x,i_y,i_z,i_m)*Bini(2,i_x,i_y,i_z,i_m)+spin(6,i_x,i_y,i_z,i_m)*Bini(3,i_x,i_y,i_z,i_m))/norm(Bini(:,i_x,i_y,i_z,i_m))
           test_ang=acos(dummy)*180/pi(1.0d0)

           if (test_ang.gt.min_angle+1.0d-4) min_angle=test_ang

           enddo
          enddo
         enddo
        enddo

       enddo


! plot the energy density
!       write(fname,'(f8.4)') kT/k_B
       toto = Trim(Adjustl(AdjustR(fname)))
       write(fname,'(a,14a,a)')'Energy_T_',(toto(i:i),i=1,len_trim(toto)),'.dat'
       write(fname2,'(a,14a,a)')'DensityOfEnergy_T_',(toto(i:i),i=1,len_trim(toto)),'.dat'

       Call EnergyDensity(fname,fname2,spin,shape_spin,tableNN,shape_tableNN,masque,shape_masque,indexNN,shape_index)

       stop 'end of the analysis of the Torques'

       end subroutine Torana

       end module m_Torana
