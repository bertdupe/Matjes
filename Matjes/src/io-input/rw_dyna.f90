      module m_dynamic
      real(kind=8) :: timestep,damping,temps,hfin(3),max_error
      integer :: duration,ti,tf,times,tstart,htimes,hstart,Efreq
      real(kind=8) :: torque_FL,torque_AFL,Ipol(3),adia,nonadia,storque,hstep(3)
      logical :: stmtorque,marche,rampe,hsweep,Ffield,i_Efield,i_torque,stmtemp
      real(kind=8), dimension(:,:,:), allocatable :: htor
! integration type
      integer :: integtype
      end module m_dynamic

      subroutine rw_dyna(N,my_lattice)
      use m_dynamic
      use m_constants
      use m_derived_types
      use m_lattice, only : spin
      implicit none
      type(lattice), intent(in) :: my_lattice
      integer, parameter :: io=23
      logical :: exists
      integer :: fin,i_x,i_y,i_z,i
      character(len=30) :: str,dummy
      integer, intent(in) :: N(3)
      real(kind=8) :: ry,xmed,ymed,ri,h,trans(3),net(3,3)

      net=my_lattice%areal

      allocate(htor(N(1),N(2),N(3)))

      i_torque=.False.
      stmtemp=.False.
      torque_FL=0.0d0
      torque_AFL=0.0d0
      Ipol=0.0d0
      adia=0.0d0
      nonadia=0.0d0
      storque=0.0d0
      hstep=0.0d0
      damping=0.0d0
      hfin=0.0d0
      htimes=0
      hstart=0
      hfin=0.0d0
      timestep=0.0d0
      damping=0.0d0
      temps=0.0d0
      stmtorque=.False.
      marche=.False.
      rampe=.False.
      hsweep=.False.
      Ffield=.False.
      i_Efield=.False.
      Efreq=1

      inquire (file='dyna.in',exist=exists)
      if (.not. exists) then
      write(6,*) 'no input file for spin dynamics'
      STOP
      endif

      open (io,file='dyna.in',form='formatted',status='old' & 
        ,action='read')

      rewind(io)
      do
      read (io,'(a)',iostat=fin) str
        if (fin /= 0) exit
        str= trim(adjustl(str))
        if (len_trim(str)==0) cycle

        if (str(1:1) == '#' ) cycle
!cccc We start to read the input
        if ( str(1:11) == 'integration') then
           backspace(io)
           read(io,*) dummy, integtype
           if (integtype.eq.5) then
            backspace(io)
            read(io,*) dummy,integtype,max_error
           endif
          endif
        if ( str(1:8) == 'timestep') then
           backspace(io)
           read(io,*) dummy, timestep
          endif
        if ( str(1:5) == 'Efreq') then
           backspace(io)
           read(io,*) dummy, Efreq
          endif
        if ( str(1:8) == 'duration') then
           backspace(io)
           read(io,*) dummy, duration
          endif
        if ( str(1:4) == 'step') then
           backspace(io)
           read(io,*) dummy, marche, ti, tf
          endif
        if ( str(1:5) == 'rampe') then
           backspace(io)
           read(io,*) dummy, rampe, ry, times, tstart
           temps=ry*k_B
          endif
        if ( str(1:7) == 'STMtemp') then
           backspace(io)
           read(io,*) dummy, stmtemp
           temps=ry*k_B
          endif
        if ( str(1:6) == 'Hsweep') then
           backspace(io)
           read(io,*) dummy, hsweep,(hstep(i),i=1,3),htimes,hstart,(hfin(i),i=1,3)
          endif
        if ( str(1:7) == 'damping') then
           backspace(io)
           read(io,*) dummy, damping
          endif
        if ( str(1:6) == 'torque') then
           backspace(io)
           read(io,*) dummy, torque_FL
           if (dabs(torque_FL).gt.1.0d-8) i_torque=.true.
          endif
        if ( str(1:10) == 'damptorque') then
           backspace(io)
           read(io,*) dummy, torque_AFL
          endif
        if ( str(1:4) == 'Ipol') then
           backspace(io)
           read(io,*) dummy, trans(:)
           Ipol=trans/(trans(1)**2+trans(2)**2+trans(3)**2)
          endif
        if ( str(1:4) == 'adia') then
           backspace(io)
           read(io,*) dummy, adia
          endif
        if ( str(1:6) == 'Ffield') then
           backspace(io)
           read(io,*) dummy, Ffield
          endif
        if ( str(1:6) == 'Efield') then
           backspace(io)
           read(io,*) dummy, i_Efield
          endif
        if ( str(1:7) == 'nonadia') then
           backspace(io)
           read(io,*) dummy, nonadia
          endif
        if ( str(1:9) == 'stmtorque') then
           backspace(io)
           read(io,*) dummy, stmtorque, storque, ri, h
          endif
      enddo 

      htor=0.0d0

      if (stmtorque) then

        write(*,*) 'SPSTM torque set up'
        xmed=dble(N(1))/2.0d0*net(1,1)+dble(N(2))/2.0d0*net(2,1)
        ymed=dble(N(1))/2.0d0*net(1,2)+dble(N(2))/2.0d0*net(2,2)

        do i_x=1,N(1)
         do i_y=1,N(2)
          do i_z=1,N(3)

        htor(i_x,i_y,i_z)=dexp(-dsqrt((xmed-spin(1,i_x,i_y,i_z,1))**2+ &
           (ymed-spin(2,i_x,i_y,i_z,1))**2+h**2)/ri)

          enddo
         enddo
        enddo
      endif
      end subroutine rw_dyna
