      module m_createspinfile
      interface CreateSpinFile
       module procedure CreateSpinFile_usernamed
       module procedure CreateSpinFile_I_simple_5d
       module procedure CreateSpinFile_R_simple_5d
       module procedure CreateSpinFile_simple_4d
       module procedure CreateSpinFile_end
       module procedure CreateSpinFile_simple_5d_usernamed
      end interface
      contains
!
! Create files to plot with Povray
!
!
! ===============================================================
      subroutine CreateSpinFile_simple_5d_usernamed(fname,signature,spin,shape_spin)
      use m_constants, only : pi
      Implicit none
      real(kind=8), intent(in) :: spin(:,:,:,:,:)
      integer, intent(in) :: shape_spin(:),signature
      character(len=*), intent(in) :: fname
!     Slope Indexes for three dim spins
      INTEGER :: i_x,i_y,i_z,i_m
! coordinate of the spins
      integer :: X,Y,Z
! position of the spins
      integer :: Rx,Ry,Rz
!     calculating the angles
      real(kind=8) :: theta, phi
!     Is the row even (1) or not (0)
      Integer :: i,phase
!     colors
      real(kind=8) :: phi_color, Delta, widthc
      real(kind=8) :: Rc,Gc,Bc
      character(len=50) :: toto,fname2,fname3

      toto=trim(adjustl(fname))
      write(fname2,'(I10)') signature
      fname3=trim(adjustl(fname2))
      write(fname2,'(50a,a,50a,a)')(toto(i:i),i=1,len_trim(toto)),'_',(fname3(i:i),i=1,len_trim(fname3)),'.dat'

      X=shape_spin(1)-3
      Y=shape_spin(1)-2
      Z=shape_spin(1)-1

      Rx=shape_spin(1)-6
      Ry=shape_spin(1)-5
      Rz=shape_spin(1)-4

!     Constants used for the color definition
      widthc=5.0d0
      Delta =PI(2.0d0/3.0d0)

      OPEN(15,FILE=fname2,status='replace')

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
          theta=acos(Spin(Z,i_x,i_y,i_z,i_m))*180.0d0/pi(1.0d0)
        else
          theta=90.0d0-dsign(90.0d0,Spin(Z,i_x,i_y,i_z,i_m))
        endif

        phi=atan2(Spin(Y,i_x,i_y,i_z,i_m),Spin(X,i_x,i_y,i_z,i_m))

        phi=phi*180.0d0/pi(1.0d0)

!       Calcualting the color as a function of the angle in or
!       out of the plane
        phi_color=pi(theta/300.0d0*2.0d0)
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

      END subroutine CreateSpinFile_simple_5d_usernamed
! ===============================================================

! ===============================================================
      subroutine CreateSpinFile_usernamed(fname,spin,shape_spin)
      use m_constants, only : pi
      Implicit none
      real(kind=8), intent(in) :: spin(:,:,:,:,:)
      integer, intent(in) :: shape_spin(:)
      character(len=*), intent(in) :: fname
!     Slope Indexes for three dim spins
      INTEGER :: i_x,i_y,i_z,i_m
! coordinate of the spins
      integer :: X,Y,Z
! position of the spins
      integer :: Rx,Ry,Rz
!     calculating the angles
      real(kind=8) :: theta, phi
!     Is the row even (1) or not (0)
      Integer :: i,phase
!     colors
      real(kind=8) :: phi_color, Delta, widthc
      real(kind=8) :: Rc,Gc,Bc
      character(len=50) :: toto

      toto=trim(adjustl(fname))

      X=shape_spin(1)-3
      Y=shape_spin(1)-2
      Z=shape_spin(1)-1

      Rx=shape_spin(1)-6
      Ry=shape_spin(1)-5
      Rz=shape_spin(1)-4

!     Constants used for the color definition
      widthc=5.0d0
      Delta =PI(2.0d0/3.0d0)

      OPEN(15,FILE=toto,status='replace')

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
          theta=acos(Spin(Z,i_x,i_y,i_z,i_m))*180.0d0/pi(1.0d0)
        else
          theta=90.0d0-dsign(90.0d0,Spin(Z,i_x,i_y,i_z,i_m))
        endif

        phi=atan2(Spin(Y,i_x,i_y,i_z,i_m),Spin(X,i_x,i_y,i_z,i_m))

        phi=phi*180.0d0/pi(1.0d0)

!       Calcualting the color as a function of the angle in or
!       out of the plane
        phi_color=pi(theta/300.0d0*2.0d0)
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

      END subroutine CreateSpinFile_usernamed
! ===============================================================

! ===============================================================
      subroutine CreateSpinFile_end(spin,shape_spin)
      use m_constants, only : pi
      Implicit none
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
      Integer :: i,phase
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
      Delta =PI(2.0d0/3.0d0)

      OPEN(15,FILE='Spinse_end.dat',status='replace')

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
          theta=acos(Spin(Z,i_x,i_y,i_z,i_m))*180.0d0/pi(1.0d0)
        else
          theta=90.0d0-dsign(90.0d0,Spin(Z,i_x,i_y,i_z,i_m))
        endif

        phi=atan2(Spin(Y,i_x,i_y,i_z,i_m),Spin(X,i_x,i_y,i_z,i_m))

        phi=phi*180.0d0/pi(1.0d0)

!       Calcualting the color as a function of the angle in or
!       out of the plane
        phi_color=pi(theta/300.0d0*2.0d0)
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
      character(len=30) :: fname,toto

     X=shape_spin(1)-3
     Y=shape_spin(1)-2
     Z=shape_spin(1)-1

!     Constants used for the color definition
      widthc=5.0d0
      Delta =PI(2.0d0/3.0d0)

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
          theta=acos(Spin(Z,i_x,i_y,i_z))*180.0d0/pi(1.0d0)
        else
          theta=90.0d0-dsign(90.0d0,Spin(Z,i_x,i_y,i_z))
        endif

        phi=atan2(Spin(Y,i_x,i_y,i_z),Spin(X,i_x,i_y,i_z))

        phi=phi*180.0d0/pi(1.0d0)

!       Calcualting the color as a function of the angle in or
!       out of the plane
        phi_color=pi(theta/300.0d0*2.0d0)
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
      use m_parameters, only : J_il
      Implicit none
      integer, intent(in) :: signature
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
      Delta =PI(2.0d0/3.0d0)

      write(fname,'(I10)') signature
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
          theta=acos(Spin(Z,i_x,i_y,i_z,i_m))*180.0d0/pi(1.0d0)
        else
          theta=90.0d0-dsign(90.0d0,Spin(Z,i_x,i_y,i_z,i_m))
        endif

        phi=atan2(Spin(Y,i_x,i_y,i_z,i_m),Spin(X,i_x,i_y,i_z,i_m))

        phi=phi*180.0d0/pi(1.0d0)

!       Calcualting the color as a function of the angle in or
!       out of the plane
        phi_color=pi(theta/300.0d0*2.0d0)
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
      END subroutine CreateSpinFile_I_simple_5d
! ===============================================================

! ===============================================================
      subroutine CreateSpinFile_R_simple_5d(signature,spin,shape_spin)
      use m_constants, only : pi
      use m_parameters, only : J_il
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
      Delta =PI(2.0d0/3.0d0)

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
          theta=acos(Spin(Z,i_x,i_y,i_z,i_m))*180.0d0/pi(1.0d0)
        else
          theta=90.0d0-dsign(90.0d0,Spin(Z,i_x,i_y,i_z,i_m))
        endif

        phi=atan2(Spin(Y,i_x,i_y,i_z,i_m),Spin(X,i_x,i_y,i_z,i_m))

        phi=phi*180.0d0/pi(1.0d0)

!       Calcualting the color as a function of the angle in or
!       out of the plane
        phi_color=pi(theta/300.0d0*2.0d0)
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
      end module m_createspinfile
