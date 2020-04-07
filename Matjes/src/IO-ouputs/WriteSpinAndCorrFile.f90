module m_write_spin
use m_get_position
use m_derived_types
use m_io_utils
use m_io_files_utils
use m_convert

  interface WriteSpinAndCorrFile
       module procedure write_usernamed
       module procedure write_general_2d
       module procedure write_general_end
       module procedure write_general_fname
       module procedure write_general_pointer_sint
       module procedure write_general_pointer_sreal
  end interface WriteSpinAndCorrFile

private
public :: WriteSpinAndCorrFile
contains

! ===============================================================
! works mainly for the GNEB
SUBROUTINE write_usernamed(spin,filename)
Implicit none
real(kind=8), intent(in) :: spin(:,:)
character(len=*) :: filename
!     Coordinates of the Spins
Integer :: io

io=open_file_write(filename)

call dump_config(io,spin)

call close_file(filename,io)

END SUBROUTINE write_usernamed
! ===============================================================

! ===============================================================
! works only with pointers of rank 1
SUBROUTINE write_general_pointer_sint(i_count,spin,filename)
use m_convert
Implicit none
type(vec_point), intent(in) :: spin(:)
integer, intent(in) :: i_count
character(len=*), intent(in) :: filename
!     Coordinates of the Spins
integer :: io
character(len=30) :: str

str=convert(filename,i_count,'.dat')
!     write Spinconfiguration in a file for STM-simulation

io=open_file_write(str)

call dump_config(io,spin)

call close_file(str,io)

END SUBROUTINE write_general_pointer_sint
! ===============================================================

! ===============================================================
! works only with pointers of rank 1
SUBROUTINE write_general_pointer_sreal(x,spin,filename)
use m_convert
Implicit none
type(vec_point), intent(in) :: spin(:)
real(kind=8), intent(in) :: x
character(len=*), intent(in) :: filename
!     Coordinates of the Spins
integer :: io
character(len=50) :: str

str=convert(filename,x,'.dat')
!     write Spinconfiguration in a file for STM-simulation

io=open_file_write(str)

call dump_config(io,spin)

call close_file(str,io)

END SUBROUTINE write_general_pointer_sreal
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
SUBROUTINE write_general_fname(fname,my_lattice)
Implicit none
type(lattice), intent(in) :: my_lattice
character(len=*), intent(in) :: fname
!     Coordinates of the Spins
integer :: io

!     write Spinconfiguration in a file for STM-simulation

io=open_file_write(fname)

call dump_config(io,my_lattice)

call close_file(fname,io)

END SUBROUTINE write_general_fname
! ===============================================================

! ===============================================================
! works only real matrix of rank 2
SUBROUTINE write_general_2d(i_count,matrix,filename)
use m_convert
Implicit none
real(kind=8), intent(in) :: matrix(:,:)
integer, intent(in) :: i_count
character(len=*), intent(in) :: filename
!     Coordinates of the Spins
integer :: io
character(len=30) :: str

str=convert(filename,i_count,'.dat')
!     write Spinconfiguration in a file for STM-simulation

io=open_file_write(str)

call dump_config(io,matrix)

call close_file(str,io)

END SUBROUTINE write_general_2d

end module m_write_spin
