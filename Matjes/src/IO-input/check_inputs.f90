      subroutine check_inputs(nlattice)
      implicit none
! intent
      integer, intent(out) :: nlattice

! local
      integer, parameter :: io=1
      logical :: exists
      character(len=100) :: str
      character(len=10) :: dummy
      integer :: fin

      nlattice=-1

! open the input
      inquire (file='input',exist=exists)
      if (.not. exists) then
      write(6,*) 'no input file'
      STOP
      endif

      open (io,file='input',form='formatted',status='old',action='read')

      rewind(io)
      do
      read (io,'(a)',iostat=fin) str
        if (fin /= 0) exit
        str= trim(adjustl(str))
        if (len_trim(str)==0) cycle

        if (str(1:1) == '#' ) cycle

        if ( str(1:8) == 'nlattice') then
           backspace(io)
           read(io,*) dummy, nlattice
        endif
      enddo
      close(io)

      end subroutine check_inputs
