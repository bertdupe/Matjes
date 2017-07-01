      subroutine rw_zdir(phase,jz,Nei_z,jil,Nei_il,dim_lat)
      implicit none
! variable shared in order to gain piece
      real (kind=8), intent(inout) :: jz(12),jil(12)
      integer, intent(inout) :: phase,Nei_z,Nei_il
      integer, intent(in) :: dim_lat(3)
! variable for the file reading
      integer, parameter :: io=24
! dumy variable
      character(len=100) :: str
      integer :: j,fin
      character(len=5) :: typ,dummy


      open (io,file='zdir.in',form='formatted',status='old' &
        ,action='read')

      rewind(io)
      do
      read (io,'(a)',iostat=fin) str
        if (fin /= 0) exit
        str= trim(adjustl(str))
        if (len_trim(str)==0) cycle

        if (str(1:1) == '#' ) cycle
! We start to read the input
        if ( str(1:5) == 'J_1sl') then
           backspace(io)
           read(io,*) dummy, jil(1)
          endif
        if ( str(1:5) == 'J_2sl') then
           backspace(io)
           read(io,*) dummy, jil(2)
          endif
        if ( str(1:5) == 'J_3sl') then
           backspace(io)
           read(io,*) dummy, jil(3)
          endif
        if ( str(1:5) == 'J_4sl') then
           backspace(io)
           read(io,*) dummy, jil(4)
          endif
        if ( str(1:5) == 'J_5sl') then
           backspace(io)
           read(io,*) dummy, jil(5)
          endif
        if ( str(1:5) == 'J_6sl') then
           backspace(io)
           read(io,*) dummy, jil(6)
          endif
        if ( str(1:5) == 'J_7sl') then
           backspace(io)
           read(io,*) dummy, jil(7)
          endif
        if ( str(1:5) == 'J_8sl') then
           backspace(io)
           read(io,*) dummy, jil(8)
          endif
        if ( str(1:5) == 'J_9sl') then
           backspace(io)
           read(io,*) dummy, jil(9)
          endif
        if ( str(1:5) == 'J_10sl') then
           backspace(io)
           read(io,*) dummy, jil(10)
          endif
        if ( str(1:5) == 'J_1z') then
           backspace(io)
           read(io,*) dummy, jz(1)
          endif
        if ( str(1:5) == 'J_2z') then
           backspace(io)
           read(io,*) dummy, jz(2)
          endif
        if ( str(1:5) == 'J_3z') then
           backspace(io)
           read(io,*) dummy, jz(3)
          endif
        if ( str(1:5) == 'J_4z') then
           backspace(io)
           read(io,*) dummy, jz(4)
          endif
        if ( str(1:5) == 'J_5z') then
           backspace(io)
           read(io,*) dummy, jz(5)
          endif
        if ( str(1:5) == 'J_6z') then
           backspace(io)
           read(io,*) dummy, jz(6)
          endif
        if ( str(1:5) == 'J_7z') then
           backspace(io)
           read(io,*) dummy, jz(7)
          endif
        if ( str(1:5) == 'J_8z') then
           backspace(io)
           read(io,*) dummy, jz(8)
          endif
        if ( str(1:5) == 'J_9z') then
           backspace(io)
           read(io,*) dummy, jz(9)
          endif
        if ( str(1:5) == 'J_10z') then
           backspace(io)
           read(io,*) dummy, jz(10)
          endif
        if ( str(1:4) == 'type') then
           backspace(io)
           read(io,*) dummy, typ
           if (typ=='super') then
            if (dim_lat(3).eq.1) then
             phase=2
            else
! the most horrible case: one more atom per unit cell and superlattice in the z direction
             phase=3
            endif
           elseif (typ=='bulk') then
            phase=1
           endif
          endif
      enddo

      Nei_il=count(abs(jil).gt.1.0d-8)
      Nei_z=count(abs(jz).gt.1.0d-8)

      end subroutine rw_zdir
