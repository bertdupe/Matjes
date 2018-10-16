module m_write_spin
use m_get_position
use m_derived_types
use m_io_utils
use m_io_files_utils
use m_convert

  interface WriteSpinAndCorrFile
       module procedure write_SD
       module procedure write_usernamed
       module procedure write_MC
       module procedure write_general_5d
       module procedure write_general_4d
       module procedure write_general_end
       module procedure write_general_pointer
  end interface WriteSpinAndCorrFile

private
public :: WriteSpinAndCorrFile
contains

! ===============================================================
! works only with matrices of rank 5 with the name given by the user
      SUBROUTINE write_usernamed(fname,spin,shape_spin)
      Implicit none
      real(kind=8), intent(in) :: spin(:,:,:,:,:)
      integer, intent(in) :: shape_spin(:)
      character(len=*) :: fname
!     Coordinates of the Spins
      Integer :: i_x,i_y,i_z,i_m,j_lat,size_z
      character(len=50) :: toto

      size_z=shape_spin(5)

!     write Spinconfiguration in a povray-file

      toto=trim(adjustl(fname))

!     write Spinconfiguration in a file for STM-simulation

      OPEN(15,FILE=toto,status='unknown')

      do i_z=1,shape_spin(4)
       do i_y=1,shape_spin(3)
        do i_x=1,shape_spin(2)

        Write(15,'(7f14.8)') ((Spin(j_lat,i_x,i_y,i_z,i_m), j_lat=1,shape_spin(1)),i_m=1,size_z)

        enddo
       enddo
      enddo

      Close(15)

      END SUBROUTINE write_usernamed
! ===============================================================

! ===============================================================
! works only with pointers of rank 1
SUBROUTINE write_general_pointer(i_count,spin)
use m_convert
Implicit none
type(vec_point), intent(in) :: spin(:)
integer, intent(in) :: i_count
!     Coordinates of the Spins
integer :: io
character(len=30) :: str

str=convert('SpinSTM_',i_count,'.dat')
!     write Spinconfiguration in a file for STM-simulation

io=open_file_write(str)

call dump_config(io,spin)

call close_file(str,io)

END SUBROUTINE write_general_pointer
! ===============================================================

! ===============================================================
! works only with matrices of rank 4
SUBROUTINE write_general_end(my_lattice)
Implicit none
type(lattice), intent(in) :: my_lattice
!     Coordinates of the Spins
integer :: io

!     write Spinconfiguration in a file for STM-simulation

io=open_file_write('SpinSTM_end.dat')

call dump_config(io,my_lattice)

call close_file('SpinSTM_end.dat',io)

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
!     Coordinates of the Spins
      Integer :: i_x,i_y,i_z,j_lat,i
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
      toto=trim(adjustl(name_in))
      fcount=trim(adjustl(fname))

!     write Spinconfiguration in a file for STM-simulation
      write(fname,'(10a,a,10a,a)')(toto(i:i),i=1,len_trim(toto)),'_',(fcount(i:i),i=1,len_trim(fcount)),'.dat'
      OPEN(15,FILE=fname)

      do i_z=1,shape_spin(4)
       do i_y=1,shape_spin(3)
        do i_x=1,shape_spin(2)

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

      do i_z=1,shape_spin(4)
       do i_y=1,shape_spin(3)
        do i_x=1,shape_spin(2)

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

      do i_z=1,shape_spin(4)
       do i_y=1,shape_spin(3)
        do i_x=1,shape_spin(2)

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

      size_z=shape_spin(5)

!     write Spinconfiguration in a povray-file
      write(fname,'(I8)') j

!     write Spinconfiguration in a file for STM-simulation

      toto=trim(adjustl(fname))
      write(fname,'(a,18a,a)')'SpinSTM_',(toto(i:i),i=1, &
     & len_trim(toto)),'.dat'
      OPEN(15,FILE=fname)

      do i_z=1,shape_spin(4)
       do i_y=1,shape_spin(3)
        do i_x=1,shape_spin(2)

         Write(15,'(7f14.8)') ((Spin(j_lat,i_x,i_y,i_z,i_m), j_lat=1,shape_spin(1)), &
     & i_m=1,size_z)

        enddo
       enddo
      enddo

      Close(15)

      END SUBROUTINE write_SD
! ===============================================================
      end module m_write_spin
