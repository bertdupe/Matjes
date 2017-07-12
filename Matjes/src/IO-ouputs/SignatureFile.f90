! ===============================================================
      SUBROUTINE SignatureFile(Inex,l1,nam,l2,pos)
      use m_parameters
      use m_efield
      use m_rw_lattice
      use m_constants
      Implicit none
!     slope variable
      integer, intent(in) :: Inex,l1,l2
      character(len=l1), intent(in) :: nam
      character(len=l2), intent(in) :: pos
      character(len=10) lattice
      integer :: i

      if (pos /= 'append') then
      OPEN(Inex,FILE=nam,action='write',status='unknown', &
       position=pos,form='formatted')
      end if
!     Vector

        Write(Inex,'(a)') '#'
        Write(Inex,'(3(a,I15,2x))') '#Nx=', dim_lat(1),' Ny=',dim_lat(2),' Nz=', dim_lat(3)
        Write(Inex,'(a,I15)') '#Total MC steps for the average (MC_Steps=)',Total_MC_Steps
        Write(Inex,'(a,I5)') '#Number of calculated Temperatures (n_TSteps=)',n_Tsteps
        Write(Inex,'(a,3f12.5)') '#External magnetic field', H_ext
        Write(Inex,'(2(a,f10.5,2x))') '#Temperature range:', kTini/k_B, '-', kTfin/k_B
!        Write(Inex,*) '#Value of cone angle:', cone
!     pause
        If (all(Periodic_log)) then
          Write(Inex,'(a)') '#with periodic boundary condition'
        else
          if (Periodic_log(1)) &
         Write(Inex,'(a)') '#with periodic boundary condition in x direction'
          if (Periodic_log(2)) &
         Write(Inex,'(a)') '#with periodic boundary condition in y direction'
          if (Periodic_log(3)) &
         Write(Inex,'(a)') '#with periodic boundary condition in z direction'
        endif
        Write(Inex,'(a,I10)') '#Nb of steps for the fetching      ', T_relax
        Write(Inex,'(a,I10)') '#Autocorrelation Time ', T_auto
        Write(Inex,'(a)') '#Exchange Interactions Jij'
        do i=2,size(J_ij(:,1)),2
!         write(Inex,1001) '#','J_',i-1,'=',J_ij(i-1)+me(i-1)*maxval(Efield_Jij),'J_',i,'=',J_ij(i)+me(i)*maxval(Efield_Jij)
         write(Inex,1001) '#','J_',i-1,'=',J_ij(i-1,1),'J_',i,'=',J_ij(i,1)
        enddo
        if (i_sliptun) then
        Write(Inex,'(a)') '# Interlayer Jijs'
         do i=2,size(J_ij(:,2)),2
         write(Inex,1001) '#','J_il',i-1,'=',J_ij(i-1,2),'J_il',i,'=',J_ij(i,2)
        enddo
        endif
        Write(Inex,'(a)') '#Exchange Interactions J_il (if necessary)'
        do i=2,size(J_il),2
         write(Inex,1001) '#','J_il',i-1,'=',J_il(i-1),'J_il',i,'=',J_il(i)
        enddo
1001    format(a,4x,2(a,I2,a,x,f12.7,2x))
        Write(Inex,'(a,3f12.7)') '#anisotropy', D_ani
        Write(Inex,'(a,3f5.2)') '#easy axis', EA
        if (i_biq) Write(Inex,'(a,f12.7)') '#J_B', J_B
        if (i_DM) Write(Inex,'(a,3f12.7)') '#DM terms', DM
        if (i_four) Write(Inex,'(a,f12.7)') '#K_1', K_1
        if (i_dip) Write(Inex,'(a)') '#dipole dipole term taken into account'
        Write(Inex,'(a,f12.7)') '#cone  ', cone
        Write(Inex,'(3(a,f6.3,3x))') '#C_jij=', c_Ji, 'C_Ki=', c_Ki, 'C_DM=', c_DM
        Write(Inex,'(2(a,f6.3,3x))') '#C_JB=', c_JB, 'C_ani=', c_ani
        if (pos /= 'append') Close(Inex)


      END SUBROUTINE SignatureFile
! ===============================================================
