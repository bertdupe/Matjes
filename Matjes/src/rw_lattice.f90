      module m_rw_lattice
      use m_derived_types
! variable that we read
      real (kind=8) :: net(3,3)
      integer :: natom
      type (cell) motif
! We are defining the world and the dimension of the problem
! the world says if you are 1, 2 or 3D
      integer, dimension(:), allocatable :: world
! N_site is the total number of site
      integer :: N_site,dim_lat(3)
! system type: n_system
! 1: 1D
! 12: 1D avec more than 1 atom per cell
! ....
      integer :: n_system
      end module m_rw_lattice

      subroutine rw_lattice()
      use m_rw_lattice
      use m_vector, only :norm
      use m_derived_types
      implicit none
! dummy variable
      integer, parameter  :: io=2
      integer :: fin,i,j,alloc_check
      logical :: exists
      character(len=100) :: str
      character(len=10) :: dummy
      real (kind=8) :: alat(3),ss

      inquire (file='lattice.in',exist=exists)
      if (.not. exists) then
      write(6,*) 'no lattice.in file'
      STOP
      endif

      open (io,file='lattice.in',form='formatted',status='old',action='read')

      rewind(io)
      do
      read (io,'(a)',iostat=fin) str
       if (fin /= 0) exit
       str= trim(adjustl(str))
       if (len_trim(str)==0) cycle

      if (str(1:1) == '#' ) cycle
! We start to read the input
       if ( str(1:4) == 'alat') then
         backspace(io)
         read(io,*) dummy,(alat(i),i=1,3)
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
           read(io,*) dummy, (dim_lat(i),i=1,3)
         if (dim_lat(1)*dim_lat(2)*dim_lat(3).eq.0) then
           write(6,*) 'dimension of lattice must be >0 along all directions'
           stop
         elseif((dim_lat(1).ne.dim_lat(2)).and.(dim_lat(2).ne.1)) then
           write(6,*) 'asymetric cell choosen'
         elseif((dim_lat(1).ne.dim_lat(2)).and.(dim_lat(2).eq.1)) then
          write(6,*) 'chain geometry choosen'
         endif
        N_site=dim_lat(1)*dim_lat(2)*dim_lat(3)
       endif
       if ( str(1:5) == 'motif') then
        backspace(io)
         read(io,*) dummy, natom, dummy
         allocate(motif%pos(natom,3),motif%mom(natom),motif%i_m(natom),motif%type(natom),stat=alloc_check)
         if (alloc_check.ne.0) write(6,'(a)') 'can not allocate motif'
          motif%pos=0.0d0
          motif%mom=0.0d0
          motif%type=1
          motif%i_m=.False.
          do i=1,natom
           read(io,*) (motif%pos(i,j),j=1,3),motif%mom(i)
          enddo
       endif
      enddo
      close(io)

      call check_motif(motif)

      do j=1,3
       net(j,:)=net(j,:)*alat(j)
      enddo

! size of the world
      if ((dim_lat(3).eq.1).and.(dim_lat(2).eq.1)) then
       allocate(world(1))
       world(1)=dim_lat(1)
       n_system=1
       if (count(motif%i_m).gt.1) n_system=12
      elseif (dim_lat(3).eq.1) then
       allocate(world(2))
       world=(/dim_lat(1),dim_lat(2)/)
       n_system=2
       if (count(motif%i_m).gt.1) n_system=22
      elseif ((dim_lat(3).eq.1).and.(dim_lat(2).eq.1).and.(dim_lat(1).eq.1)) then
       write(6,*) "dimension of the problem not correct"
       stop
      else
       allocate(world(3))
       world=(/dim_lat(1),dim_lat(2),dim_lat(3)/)
       n_system=3
       if (count(motif%i_m).gt.1) n_system=32
      endif
! takes into account that we might have more than one magnetic atom per unit cell

      call update_lattice(motif,N_site)

      end subroutine rw_lattice
