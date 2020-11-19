module m_createspinfile
use m_constants, only : pi
use m_get_position
use m_derived_types
use m_io_files_utils
use m_io_utils
use m_convert

  interface CreateSpinFile
       module procedure CreateSpinFile_lattice_usernamed
       module procedure CreateSpinFile_I_simple_5d
       module procedure CreateSpinFile_I_lattice_sint
       module procedure CreateSpinFile_I_lattice_sreal
       module procedure CreateSpinFile_R_simple_5d
       module procedure CreateSpinFile_simple_4d
       module procedure CreateSpinFile_end
       module procedure CreateSpinFile_2d_sint
       module procedure CreateSpinFile_2d_sreal
       module procedure CreateSpinFile_usernamed_spin
  end interface

private
public :: CreateSpinFile

contains
!
! Create files to plot with Povray
!
!

! ===============================================================
subroutine CreateSpinFile_2d_sint(fname,signature,spin)
Implicit none
real(kind=8), intent(in) :: spin(:,:)
integer, intent(in) :: signature
character(len=*), intent(in) :: fname
!     Slope Indexes for three dim spins
character(len=50) :: name

name=convert(fname,signature,'.dat')

call CreateSpinFile_usernamed_spin(spin,name)

END subroutine CreateSpinFile_2d_sint
! ===============================================================

! ===============================================================
subroutine CreateSpinFile_2d_sreal(fname,signature,spin)
Implicit none
real(kind=8), intent(in) :: spin(:,:)
real(kind=8), intent(in) :: signature
character(len=*), intent(in) :: fname
!     Slope Indexes for three dim spins
character(len=50) :: name

name=convert(fname,signature,'.dat')

call CreateSpinFile_usernamed_spin(spin,name)

END subroutine CreateSpinFile_2d_sreal
! ===============================================================

! ===============================================================
subroutine CreateSpinFile_usernamed_spin(spin,fname)
Implicit none
character(len=*), intent(in) :: fname
real(kind=8), intent(in) :: spin(:,:)
! coordinate of the spins
integer :: io

io=open_file_write(trim(adjustl(fname)))

call dump_spinse(io,spin)

call close_file(trim(adjustl(fname)),io)

END subroutine CreateSpinFile_usernamed_spin
! ===============================================================

! ===============================================================
subroutine CreateSpinFile_lattice_usernamed(fname,my_lattice,my_motif)
Implicit none
character(len=*), intent(in) :: fname
type(lattice), intent(in) :: my_lattice
type(t_cell), intent(in) :: my_motif
! lattice of the positions
real(kind=8), allocatable, dimension(:,:,:,:,:) :: position
! coordinate of the spins
integer :: N(4),Natom_motif,io
! calculating the angles
real(kind=8) :: r(3,3)


N=shape(my_lattice%ordpar%l_modes)
Natom_motif=count(my_motif%atomic(:)%moment.gt.0.0d0)
r=my_lattice%areal

allocate(position(3,N(1),N(2),N(3),N(4)))
position=0.0d0

call get_position(position,N,r,my_motif)

io=open_file_write(trim(adjustl(fname)))

call dump_spinse(io,my_lattice,position)

call close_file(trim(adjustl(fname)),io)

deallocate(position)

END subroutine CreateSpinFile_lattice_usernamed
! ===============================================================

! ===============================================================
subroutine CreateSpinFile_end(my_lattice,my_motif)
Implicit none
type(lattice), intent(in) :: my_lattice
type(t_cell), intent(in) :: my_motif
! lattice of the positions
real(kind=8), allocatable, dimension(:,:,:,:,:) :: position
! coordinate of the spins
integer :: N(4),Natom_motif,io
! calculating the angles
real(kind=8) :: r(3,3)

N=shape(my_lattice%ordpar%l_modes)
Natom_motif=count(my_motif%atomic(:)%moment.gt.0.0d0)
r=my_lattice%areal

allocate(position(3,N(1),N(2),N(3),N(4)))
position=0.0d0

call get_position(position,N,r,my_motif)

! dump the positions of the sites

io=open_file_write('positions.dat')

call dump_config(io,position)

call close_file('positions.dat',io)

! dump the povray file

io=open_file_write('Spinse_end.dat')

call dump_spinse(io,my_lattice,position)

call close_file('Spinse_end.dat',io)

deallocate(position)

END subroutine CreateSpinFile_end
! ===============================================================

! ===============================================================
      subroutine CreateSpinFile_simple_4d(name_in,spin,shape_spin)
      use m_constants, only : pi
      Implicit none
      character(len=*), intent(in) :: name_in
      real(kind=8), intent(in) :: spin(:,:,:,:)
      integer, intent(in) :: shape_spin(:)
!     Slope Indexes for three dim spins
      INTEGER :: i_x,i_y,i_z
! coordinate of the spins
      integer :: X,Y,Z

!     calculating the angles
      real(kind=8) :: theta, phi
!     Is the row even (1) or not (0)
      Integer :: i
!     colors
      real(kind=8) :: phi_color, Delta, widthc
      real(kind=8) :: Rc,Gc,Bc
!   name of files
      character(len=50) :: fname,toto

     X=shape_spin(1)-3
     Y=shape_spin(1)-2
     Z=shape_spin(1)-1

!     Constants used for the color definition
      widthc=5.0d0
      Delta =PI*2.0d0/3.0d0

      toto=trim(adjustl(name_in))
      write(fname,'(a,18a,a)')'Spinse_',(toto(i:i),i=1,len_trim(toto)),'.dat'
      OPEN(15,FILE=fname,status='unknown')

      Do i_z=1,shape_spin(4)
       Do i_y=1,shape_spin(3)
        Do i_x=1,shape_spin(2)

!       Yes for these formulars it is helpfull to make a picture.
!       The initial object is a cone with its top at
!       coordinates (0,0,1). First turn it around the y-achse into
!       the x-z-plane around angly, then turn it around the z-achse
!       into the right position.
!       Then translate it to the right r_x,r_y position.
        if (abs(Spin(Z,i_x,i_y,i_z)).lt.1.0d0) then
          theta=acos(Spin(Z,i_x,i_y,i_z))*180.0d0*pi
        else
          theta=90.0d0-dsign(90.0d0,Spin(Z,i_x,i_y,i_z))
        endif

        phi=atan2(Spin(Y,i_x,i_y,i_z),Spin(X,i_x,i_y,i_z))

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

         write(15,'(2(a,f16.8),3(a,I8),3(a,f16.8),a)') 'Spin(', &
     & theta,',',phi,',',i_x,',',i_y,',',i_z,',', &
     & Rc,',',Bc,',',Gc,')'

         enddo
        enddo
       enddo

      Close(15)
      END subroutine CreateSpinFile_simple_4d
! ===============================================================

! ===============================================================
subroutine CreateSpinFile_I_simple_5d(signature,spin,shape_spin)
use m_constants, only : pi
Implicit none
integer, intent(in) :: signature,shape_spin(:)
real(kind=8), intent(in) :: spin(:,:,:,:,:)
! internal variables
!     Slope Indexes for three dim spins
INTEGER :: i_x,i_y,i_z,i_m
!     calculating the angles
real(kind=8) :: theta, phi
!     Is the row even (1) or not (0)
Integer :: i
!     colors
real(kind=8) :: Rc,Gc,Bc
! position
real(kind=8) :: pos(3)
!   name of files
character(len=50) :: fname,toto

write(fname,'(I10)') signature
toto=trim(adjustl(fname))
write(fname,'(a,18a,a)')'Spinse_',(toto(i:i),i=1,len_trim(toto)),'.dat'
OPEN(15,FILE=fname,status='unknown')

do i_m=1,shape_spin(5)
   Do i_z=1,shape_spin(4)
      Do i_y=1,shape_spin(3)
         Do i_x=1,shape_spin(2)


         call get_colors(Rc,Gc,Bc,theta,phi,spin(4:6,i_x,i_y,i_z,i_m))

        pos=spin(1:3,i_x,i_y,i_z,i_m)

         write(15,'(8(a,f16.8),a)') 'Spin(', &
     & theta,',',phi,',',pos(1),',',pos(2),',',pos(3),',',Rc,',',Bc,',',Gc,')'

         enddo
      enddo
   enddo
enddo

Close(15)
END subroutine CreateSpinFile_I_simple_5d
! ===============================================================

! ===============================================================
subroutine CreateSpinFile_I_lattice_sint(signature,spin)
Implicit none
integer, intent(in) :: signature
type(vec_point), intent(in) :: spin(:)
! internal variables
Integer :: io
!   name of files
character(len=50) :: fname

fname=convert('Spinse_',signature,'.dat')
io=open_file_write(fname)

call dump_spinse(io,spin)

call close_file(fname,io)

END subroutine CreateSpinFile_I_lattice_sint
! ===============================================================

! ===============================================================
subroutine CreateSpinFile_I_lattice_sreal(signature,spin)
Implicit none
real(kind=8), intent(in) :: signature
type(vec_point), intent(in) :: spin(:)
! internal variables
Integer :: io
!   name of files
character(len=50) :: fname

fname=convert('Spinse_',signature,'.dat')
io=open_file_write(fname)

call dump_spinse(io,spin)

call close_file(fname,io)

END subroutine CreateSpinFile_I_lattice_sreal
! ===============================================================

! ===============================================================
      subroutine CreateSpinFile_R_simple_5d(signature,spin,shape_spin)
      use m_constants, only : pi
      Implicit none
      real(kind=8), intent(in) :: signature
      real(kind=8), intent(in) :: spin(:,:,:,:,:)
      integer, intent(in) :: shape_spin(:)
!     Slope Indexes for three dim spins
      INTEGER :: i_x,i_y,i_z,i_m
! coordinate of the spins
      integer :: X,Y,Z
! position of the spins
      integer :: Rx,Ry,Rz
!     calculating the angles
      real(kind=8) :: theta, phi
!     Is the row even (1) or not (0)
      Integer :: i
!     colors
      real(kind=8) :: phi_color, Delta, widthc
      real(kind=8) :: Rc,Gc,Bc
!   name of files
      character(len=30) :: fname,toto

      X=shape_spin(1)-3
      Y=shape_spin(1)-2
      Z=shape_spin(1)-1

      Rx=shape_spin(1)-6
      Ry=shape_spin(1)-5
      Rz=shape_spin(1)-4

!     Constants used for the color definition
      widthc=5.0d0
      Delta =PI*2.0d0/3.0d0

      write(fname,'(f8.4)') signature
      toto=trim(adjustl(fname))
      write(fname,'(a,18a,a)')'Spinse_',(toto(i:i),i=1,len_trim(toto)),'.dat'
      OPEN(15,FILE=fname,status='unknown')

      do i_m=1,shape_spin(5)
       Do i_z=1,shape_spin(4)
        Do i_y=1,shape_spin(3)
         Do i_x=1,shape_spin(2)

!       Yes for these formulars it is helpfull to make a picture.
!       The initial object is a cone with its top at
!       coordinates (0,0,1). First turn it around the y-achse into
!       the x-z-plane around angly, then turn it around the z-achse
!       into the right position.
!       Then translate it to the right r_x,r_y position.
        if (abs(Spin(Z,i_x,i_y,i_z,i_m)).lt.1.0d0) then
          theta=acos(Spin(Z,i_x,i_y,i_z,i_m))*180.0d0/pi
        else
          theta=90.0d0-dsign(90.0d0,Spin(Z,i_x,i_y,i_z,i_m))
        endif

        phi=atan2(Spin(Y,i_x,i_y,i_z,i_m),Spin(X,i_x,i_y,i_z,i_m))

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

         write(15,'(8(a,f16.8),a)') 'Spin(', &
     & theta,',',phi,',',Spin(Rx,i_x,i_y,i_z,i_m),',',Spin(Ry,i_x,i_y,i_z,i_m),',',Spin(Rz,i_x,i_y,i_z,i_m),',', &
     & Rc,',',Bc,',',Gc,')'

         enddo
        enddo
       enddo
      enddo

      Close(15)
      END subroutine CreateSpinFile_R_simple_5d

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

end module m_createspinfile
