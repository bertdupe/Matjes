      module m_efield
      real(kind=8), allocatable :: Efield_Jij(:,:,:)
      real(kind=8) :: me(12)
      end module m_efield

      subroutine rw_efield(N,net)
      use m_efield
      use m_lattice, only : spin
      use m_vector, only : norm
      implicit none
! entry variable
      integer, intent(in) :: N(3)
      real(kind=8), intent(in) :: net(3,3)
! internal variable
      integer, parameter  :: io=1
      character(len=30) :: str
      character(len=3) :: tag,dummy
      real(kind=8) :: EF,radius,height,xmed,ymed
      integer :: k,i,fin,j
      real(kind=8) :: rx,ry,rz,rev
      logical :: exists

      allocate (Efield_Jij(N(1),N(2),N(3)))

      me=0.0d0
      Efield_Jij=0.0d0

      inquire (file='efield.in',exist=exists)
      if (exists) then
      write(6,*) 'electric field set up'

      open (io,file='efield.in',form='formatted',status='old', &
       action='read')
      rewind(io)

      do
      read (io,'(a)',iostat=fin) str
        if (fin /= 0) exit
        str= trim(adjustl(str))
        if (len_trim(str)==0) cycle

        if (str(1:1) == '#' ) cycle
! shape of the electric field
      if (str(1:5) == "shape") then
       backspace(io)
       read(io,*) dummy, tag
! if circ then read radius of the cylinder
! the field is then equal to E in the cylinder and 0 outside
       if (tag == "cir") then
        backspace(io)
        read(io,*) dummy, tag, radius, xmed, ymed
! if tip then read height
       elseif (tag == "tip") then
        backspace(io)
        read(io,*) dummy, tag, height, xmed, ymed
       endif
      xmed=norm(xmed*net(1,:))
      ymed=norm(ymed*net(2,:))
      endif
!read the field in V/A
      if (str(1:5) == "field") then
       backspace(io)
       read(io,*) dummy, EF
      endif
! read de a and b coeff
      if (str(1:5) == 'me_1') then
       backspace(io)
       read(io,*) dummy, me(1)
      endif
      if (str(1:5) == 'me_2') then
       backspace(io)
       read(io,*) dummy, me(2)
      endif
      if (str(1:5) == 'me_3') then
       backspace(io)
       read(io,*) dummy, me(3)
      endif
      if (str(1:5) == 'me_4') then
       backspace(io)
       read(io,*) dummy, me(4)
      endif
      if (str(1:5) == 'me_5') then
       backspace(io)
       read(io,*) dummy, me(5)
      endif
      if (str(1:5) == 'me_6') then
       backspace(io)
       read(io,*) dummy, me(6)
      endif
      if (str(1:5) == 'me_7') then
       backspace(io)
       read(io,*) dummy, me(7)
      endif
      if (str(1:5) == 'me_8') then
       backspace(io)
       read(io,*) dummy, me(8)
      endif
      if (str(1:5) == 'me_9') then
       backspace(io)
       read(io,*) dummy, me(9)
      endif
      if (str(1:5) == 'me_10') then
       backspace(io)
       read(io,*) dummy, me(10)
      endif
      if (str(1:5) == 'me_11') then
       backspace(io)
       read(io,*) dummy, me(11)
      endif
      if (str(1:5) == 'me_12') then
       backspace(io)
       read(io,*) dummy, me(12)
      endif
      end do
      close(io)

!!! the plane means uniforma field
       if (tag == "pla") then
        Efield_Jij=EF       
!! tip in a very simple model
       elseif (tag == "tip") then

        do i=1,N(1)
         do j=1,N(2)
          do k=1,N(3)

        Efield_Jij(i,j,k)=EF*height/(height**2+((xmed-Spin(1,i,j,k,1))**2+(ymed-Spin(2,i,j,k,1))**2))**2
          enddo
         enddo
        enddo
!! cylinder in an even simpler model
       elseif (tag == "cir") then

        do i=1,N(1)
         do j=1,N(2)
          do k=1,N(3)
          if (dsqrt((xmed-Spin(1,i,j,k,1))**2+(ymed-Spin(2,i,j,k,1))**2).lt.radius) then
           Efield_Jij(i,j,k)=EF
           else
           Efield_Jij(i,j,k)=0.0d0
          endif
          enddo
         enddo
        enddo
       endif
! if no field then
      else
       Efield_Jij=0.0d0
       me=0.0d0
      endif

      end subroutine rw_efield
