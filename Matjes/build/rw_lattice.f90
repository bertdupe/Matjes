      subroutine rw_lattice(my_lattice)
      use m_vector, only : norm,cross
      use m_constants, only : pi
      use m_derived_types
      implicit none
! out
      type(lattice), intent(out) :: my_lattice
! dummy variable
      integer, parameter  :: io=2
      integer :: fin,i,j,dimension(3)
      integer :: N_site
      logical :: Periodic_log(3)
      character(len=100) :: str
      character(len=10) :: dummy
      real (kind=8) :: a(3),ss,net(3,3),kv0(3,3)

      open (io,file='input',form='formatted',status='old',action='read')

      rewind(io)
      do
      read (io,'(a)',iostat=fin) str
       if (fin /= 0) exit
       str= trim(adjustl(str))
       if (len_trim(str)==0) cycle

      if (str(1:1) == '#' ) cycle
! We start to read the input
       if ( str(1:12) == 'Periodic_log') then
         backspace(io)
         read(io,*) dummy, Periodic_log(1), Periodic_log(2), Periodic_log(3)
       endif
       if ( str(1:4) == 'alat') then
         backspace(io)
         read(io,*) dummy,(a(i),i=1,3)
       endif
       if ( str(1:7) == 'lattice') then
         do j=1,3
          read(io,*) (net(j,i),i=1,3)
         enddo
         do j=1,3
          ss=norm(net(j,:))
          if (ss.lt.1.0d-8) stop 'lattice vectors are crap'
          net(j,:)=net(j,:)/ss
         enddo
       endif
       if ( str(1:5) == 'Nsize') then
         backspace(io)
           read(io,*) dummy, (dimension(i),i=1,3)
         if (dimension(1)*dimension(2)*dimension(3).eq.0) then
           write(6,*) 'dimension of lattice must be >0 along all directions'
           stop
         elseif((dimension(1).ne.dimension(2)).and.(dimension(2).ne.1)) then
           write(6,*) 'asymetric cell choosen'
         elseif((dimension(1).ne.dimension(2)).and.(dimension(2).eq.1)) then
          write(6,*) 'chain geometry choosen'
         endif
        N_site=dimension(1)*dimension(2)*dimension(3)
       endif

      enddo
      close(io)

      do j=1,3
       net(j,:)=net(j,:)*a(j)
      enddo

!! put the magnetic lattice in order
      my_lattice%dim_lat=dimension
      my_lattice%areal=net
      my_lattice%alat=a
      my_lattice%boundary=Periodic_log

! build up the reciprocal lattice vectors

      kv0(1,:) = pi(2.0d0)*cross(net(2,:),net(3,:))/dot_product(net(1,:),cross(net(2,:),net(3,:)))
      kv0(2,:) = pi(2.0d0)*cross(net(3,:),net(1,:))/dot_product(net(1,:),cross(net(2,:),net(3,:)))
      kv0(3,:) = pi(2.0d0)*cross(net(1,:),net(2,:))/dot_product(net(1,:),cross(net(2,:),net(3,:)))

      my_lattice%astar=kv0

      end subroutine rw_lattice
