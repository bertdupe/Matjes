      module m_write_spin
      interface WriteSpinAndCorrFile
       module procedure write_SD
       module procedure write_usernamed
       module procedure write_MC
       module procedure write_general_5d
       module procedure write_general_4d
       module procedure write_general_end
      end interface WriteSpinAndCorrFile

      contains

! ===============================================================
! works only with matrices of rank 5 with the name given by the user
      SUBROUTINE write_usernamed(fname,spin,shape_spin)
      Implicit none
      real(kind=8), intent(in) :: spin(:,:,:,:,:)
      integer, intent(in) :: shape_spin(:)
      character(len=*) :: fname
!     slope index for corre
      Integer :: i
!     Coordinates of the Spins
      Integer :: i_x,i_y,i_z,i_m,j_lat,size_z
      character(len=50) :: toto

      size_z=shape_spin(5)

!     write Spinconfiguration in a povray-file

      toto=trim(adjustl(fname))

!     write Spinconfiguration in a file for STM-simulation

      OPEN(15,FILE=toto,status='unknown')

      do i_x=1,shape_spin(2)
       do i_y=1,shape_spin(3)
        do i_z=1,shape_spin(4)

        Write(15,'(7f14.8)') ((Spin(j_lat,i_x,i_y,i_z,i_m), j_lat=1,shape_spin(1)),i_m=1,size_z)

        enddo
       enddo
      enddo

      Close(15)

      END SUBROUTINE write_usernamed
! ===============================================================

! ===============================================================
! works only with matrices of rank 4
      SUBROUTINE write_general_end(spin,shape_spin)
      Implicit none
      real(kind=8), intent(in) :: spin(:,:,:,:,:)
      integer, intent(in) :: shape_spin(:)
!     slope index for corre
      Integer :: i
!     Coordinates of the Spins
      Integer :: i_x,i_y,i_z,i_m,j_lat,size_z

      size_z=shape_spin(5)

!     write Spinconfiguration in a file for STM-simulation

      OPEN(15,FILE='SpinSTM_end.dat')

      do i_x=1,shape_spin(2)
       do i_y=1,shape_spin(3)
        do i_z=1,shape_spin(4)

        Write(15,'(7f14.8)') ((Spin(j_lat,i_x,i_y,i_z,i_m), j_lat=1,shape_spin(1)),i_m=1,size_z)

        enddo
       enddo
      enddo

      Close(15)

      END SUBROUTINE write_general_end
! ===============================================================

! ===============================================================
! works only with matrices of rank 4
      SUBROUTINE write_general_4d(name_in,i_count,spin,shape_spin)
      Implicit none
! name of the ouput
      character(LEN=*), intent(in) :: name_in
!      integer, intent(in) :: n_Tsteps
      real(kind=8), intent(in) :: spin(:,:,:,:)
! input counter for the numbering of the files
      integer, intent(in) :: i_count,shape_spin(:)
!     slope index for corre
      Integer :: i
!     Helping Variables for Modulus
      Integer :: dist_gra, add_cor, add_gra
!     Coordinates of the Spins
      Integer :: i_x,i_y,i_z,j_lat,size_z
!     name of the file
      character(len=30) :: fname,toto,fcount

!     In dependence on the number of files we
!     calculate, how many temperature steps lies between
!     to files
!     We want badly to have a picture of Spins and Correlation
!     at the lowest temperature. Nothing else is that important
!     Therefore we need that add-values
!      add_gra=mod(n_Tsteps,n_spingra)
!      add_cor=mod(n_Tsteps,n_cor)

!     write Spinconfiguration in a povray-file
      write(fname,'(I10)') i_count

!     write Spinconfiguration in a file for STM-simulation
      write(fname,'(10a,a,10a,a)')(toto(i:i),i=1,len_trim(toto)),'_',(fcount(i:i),i=1,len_trim(fcount)),'.dat'
      OPEN(15,FILE=fname)

      do i_x=1,shape_spin(2)
       do i_y=1,shape_spin(3)
        do i_z=1,shape_spin(4)

        Write(15,'(7f14.8)') (Spin(j_lat,i_x,i_y,i_z), j_lat=1,shape_spin(1))

        enddo
       enddo
      enddo

      Close(15)

      END SUBROUTINE write_general_4d
! ===============================================================

! ===============================================================
! works only with matrices of rank 5
      SUBROUTINE write_general_5d(name_in,i_count,spin,shape_spin)
      Implicit none
! name of the ouput
      character(len=*), intent(in) :: name_in
!      integer, intent(in) :: n_Tsteps
      real(kind=8), intent(in) :: spin(:,:,:,:,:)
! input counter for the numbering of the files
      integer, intent(in) :: i_count,shape_spin(:)
!     slope index for corre
      Integer :: i
!     Helping Variables for Modulus
      Integer :: dist_gra, add_cor, add_gra
!     Coordinates of the Spins
      Integer :: i_x,i_y,i_z,i_m,j_lat,size_z
!     name of the file
      character(len=30) :: fname,toto,fcount

!     In dependence on the number of files we
!     calculate, how many temperature steps lies between
!     to files
!     We want badly to have a picture of Spins and Correlation
!     at the lowest temperature. Nothing else is that important
!     Therefore we need that add-values
!      add_gra=mod(n_Tsteps,n_spingra)
!      add_cor=mod(n_Tsteps,n_cor)

      size_z=shape_spin(5)

!     write Spinconfiguration in a povray-file
      write(fname,'(I10)') i_count

!     write Spinconfiguration in a file for STM-simulation
      toto=trim(adjustl(name_in))
      fcount=trim(adjustl(fname))
      write(fname,'(10a,a,10a,a)')(toto(i:i),i=1,len_trim(toto)),'_',(fcount(i:i),i=1,len_trim(fcount)),'.dat'
      OPEN(15,FILE=fname)

      do i_x=1,shape_spin(2)
       do i_y=1,shape_spin(3)
        do i_z=1,shape_spin(4)

        Write(15,'(7f14.8)') ((Spin(j_lat,i_x,i_y,i_z,i_m), j_lat=1,shape_spin(1)), &
     & i_m=1,size_z)

        enddo
       enddo
      enddo

      Close(15)

      END SUBROUTINE write_general_5d
! ===============================================================


! ===============================================================
      SUBROUTINE write_MC(kt,spin,shape_spin)
      use m_constants, only : k_B
      Implicit none
!      integer, intent(in) :: n_kT,n_Tsteps
      real(kind=8), intent(in) :: kt,spin(:,:,:,:,:)
      integer, intent(in) :: shape_spin(:)
!     slope index for corre
      Integer :: i
!     Helping Variables for Modulus
      Integer :: dist_gra, add_cor, add_gra
!     Coordinates of the Spins
      Integer :: i_x,i_y,i_z,i_m,j_lat,size_z
!     name of the file
      character(len=30) :: fname,toto

!     In dependence on the number of files we
!     calculate, how many temperature steps lies between
!     to files
!     We want badly to have a picture of Spins and Correlation
!     at the lowest temperature. Nothing else is that important
!     Therefore we need that add-values
!      add_gra=mod(n_Tsteps,n_spingra)
!      add_cor=mod(n_Tsteps,n_cor)

      size_z=shape_spin(5)

!     write Spinconfiguration in a file for STM-simulation
      write(fname,'(f8.4)') kT/k_B
      toto=trim(adjustl(fname))
      write(fname,'(a,18a,a)')'SpinSTM_T',(toto(i:i),i=1, &
     & len_trim(toto)),'.dat'
      OPEN(15,FILE=fname)

      do i_x=1,shape_spin(2)
       do i_y=1,shape_spin(3)
        do i_z=1,shape_spin(4)

         Write(15,'(7f14.8)') ((Spin(j_lat,i_x,i_y,i_z,i_m), j_lat=1,shape_spin(1)), &
     & i_m=1,size_z)

        enddo
       enddo
      enddo

      Close(15)


      END SUBROUTINE write_MC
! ===============================================================

! ===============================================================
      SUBROUTINE write_SD(j,spin,shape_spin)
      Implicit none
!      integer, intent(in) :: n_kT,n_Tsteps
      real(kind=8), intent(in) :: spin(:,:,:,:,:)
      integer, intent(in) :: j,shape_spin(:)
!     slope index for corre
      Integer :: i
!     Helping Variables for Modulus
      Integer :: dist_gra, add_cor, add_gra
!     Coordinates of the Spins
      Integer :: i_x,i_y,i_z,i_m,j_lat,size_z
!     name of the file
      character(len=30) :: fname,toto

!     In dependence on the number of files we
!     calculate, how many temperature steps lies between
!     to files
!     We want badly to have a picture of Spins and Correlation
!     at the lowest temperature. Nothing else is that important
!     Therefore we need that add-values
!      add_gra=mod(n_Tsteps,n_spingra)
!      add_cor=mod(n_Tsteps,n_cor)

      size_z=shape_spin(5)

!     write Spinconfiguration in a povray-file
      write(fname,'(I8)') j

!     write Spinconfiguration in a file for STM-simulation

      toto=trim(adjustl(fname))
      write(fname,'(a,18a,a)')'SpinSTM_',(toto(i:i),i=1, &
     & len_trim(toto)),'.dat'
      OPEN(15,FILE=fname)

      do i_x=1,shape_spin(2)
       do i_y=1,shape_spin(3)
        do i_z=1,shape_spin(4)

         Write(15,'(7f14.8)') ((Spin(j_lat,i_x,i_y,i_z,i_m), j_lat=1,shape_spin(1)), &
     & i_m=1,size_z)

        enddo
       enddo
      enddo

      Close(15)

      END SUBROUTINE write_SD
! ===============================================================
      end module m_write_spin
