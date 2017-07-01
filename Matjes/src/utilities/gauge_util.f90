      module m_gauge
      double precision :: VectPolarRef(3)
      end module

      subroutine rw_gauge()
      use m_gauge
      implicit none
! internal variable
      real(kind=8) :: VectZ(3)
      logical :: exists
      integer, parameter  :: io=1
      integer :: fin,i
      character(len=20) :: str,dummy

      VectZ=0.0d0

      inquire (file='gauge.in',exist=exists)
      if (.not. exists) then
      write(6,*) 'no gauge.in file'
      STOP
      endif

      open (io,file='gauge.in',form='formatted',status='old', &
        action='read')

      rewind(io)
      do
      read (io,'(a)',iostat=fin) str
        if (fin /= 0) exit
        str= trim(adjustl(str))
        if (len_trim(str)==0) cycle

        if (str(1:1) == '#' ) cycle
!cccc We start to read the input
        if ( str(1:5) == 'VectZ') then
           backspace(io)
           read(io,*) dummy, (VectZ(i),i=1,3)
          endif
      enddo

      VectPolarRef=VectZ

      end subroutine rw_gauge
