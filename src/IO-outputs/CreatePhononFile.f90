module m_createphononfile
use m_get_position
use m_derived_types
use m_io_files_utils
use m_io_utils
use m_convert
implicit none

  interface CreatePhononFile
       module procedure CreatePhononFile_sig_orderpar
       module procedure CreatePhononFile_2d_sint
       module procedure CreatePhononFile_2d_sreal
       module procedure CreatePhononFile_usernamed
       module procedure CreatePhononFile_orderpar_name
  end interface

private
public :: CreatePhononFile

contains
!
! Create files to plot with Povray
!
!

! ===============================================================


subroutine CreatePhononFile_sig_orderpar(signature,ordpar)
    type(order_par), intent(in) :: ordpar
    integer,intent(in)          :: signature
    character(len=50) :: fname

    fname=convert('Povphonon_',signature,'.dat')
    Call CreatePhononFile_usernamed(ordpar%modes_v,fname)
end subroutine

subroutine CreatePhononFile_orderpar_name(ordpar,fname)
    type(order_par), intent(in) :: ordpar
    character(len=*), intent(in) :: fname
    Call CreatePhononFile_usernamed(ordpar%modes_v,fname)
END subroutine


subroutine CreatePhononFile_2d_sint(fname,signature,phonon)
Implicit none
real(kind=8), intent(in) :: phonon(:,:)
integer, intent(in) :: signature
character(len=*), intent(in) :: fname
!     Slope Indexes for three dim spins
character(len=50) :: name

name=convert(fname,signature,'.dat')

call CreatePhononFile_usernamed(phonon,name)

END subroutine CreatePhononFile_2d_sint
! ===============================================================

! ===============================================================
subroutine CreatePhononFile_2d_sreal(fname,signature,phonon)
Implicit none
real(kind=8), intent(in) :: phonon(:,:)
real(kind=8), intent(in) :: signature
character(len=*), intent(in) :: fname
!     Slope Indexes for three dim spins
character(len=50) :: name

name=convert(fname,signature,'.dat')

call CreatePhononFile_usernamed(phonon,name)

END subroutine CreatePhononFile_2d_sreal
! ===============================================================

! ===============================================================
subroutine CreatePhononFile_usernamed(phonon,fname)
Implicit none
character(len=*), intent(in) :: fname
real(kind=8), intent(in) :: phonon(:,:)
! coordinate of the spins
integer :: io

io=open_file_write(trim(adjustl(fname)))

call dump_spinse(io,phonon)

call close_file(trim(adjustl(fname)),io)

END subroutine CreatePhononFile_usernamed
! ===============================================================


! ===============================================================
! function that gets the RGB colors
! ===============================================================
subroutine get_colors(Rc,Gc,Bc,theta,phi,mode)
use m_constants,only : pi
implicit none
real(kind=8),intent(out) :: Rc,Gc,Bc,theta,phi
real(kind=8),intent(in) :: mode(3)
! internal
real(kind=8) :: widthc,Delta,phi_color

widthc=5.0d0
Delta =PI*2.0d0/3.0d0

!       Yes for these formulars it is helpfull to make a picture.
!       The initial object is a cone with its top at
!       coordinates (0,0,1). First turn it around the y-achse into
!       the x-z-plane around angly, then turn it around the z-achse
!       into the right position.
!       Then translate it to the right r_x,r_y position.

if (abs(mode(3)).lt.1.0d0) then
  theta=acos(mode(3))*180.0d0/pi
else
  theta=90.0d0-dsign(90.0d0,mode(3))
endif

phi=atan2(mode(2),mode(1))

phi=phi*180.0d0/pi

!       Calcualting the color as a function of the angle in or
!       out of the plane
phi_color=pi*theta/300.0d0*2.0d0
Rc = widthc*(cos(phi_color+0*Delta))
if (Rc.lt.0.000001d0)  Rc=0.0d0
Gc = widthc*(cos(phi_color+1*Delta))
if (Gc.lt.0.000001d0)  Gc=0.0d0
Bc = widthc*(cos(phi_color+2*Delta))
if (Bc.lt.0.000001d0)  Bc=0.0d0

end subroutine get_colors
! ===============================================================

end module m_createphononfile
