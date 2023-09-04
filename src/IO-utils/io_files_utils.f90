module m_io_files_utils
use,intrinsic :: iso_fortran_env, only : output_unit, error_unit
implicit none

private
public :: open_file_read,close_file,open_file_write

contains

!!!!!!!!!!!!!!!!!!!!!!!!
! part for pure I/O
!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!
! PUBLIC
!!!!!!!!!!!!!!!!!!!!!!!!

! function that open a file and return the unit for reading
integer function open_file_write(fname,stats,accesss)
implicit none
character(len=*), intent(in) :: fname
character(len=*), intent(in), optional :: stats,accesss
! internal variables
logical :: i_file
integer :: io,error
character(len=30) :: stat_int,access_int

stat_int='unknown'
access_int='direct'
io=-1

if (present(stats)) then
  stat_int=trim(stats)
  select case (trim(stats))
    case('unknown')
      i_file=inquire_file(fname)
    case('old')
      i_file=.true.
    case('new')
      i_file=.true.
    case default
      i_file=inquire_file(fname)
  end select
else
  stat_int='unknown'
  i_file=inquire_file(fname)
endif

if (present(accesss)) access_int=trim(accesss)

if (i_file) then
! find unit of the file that is alread
  io=inquire_file_open(fname)
endif

if (io.lt.0) io=inquire_file_unit(io)

open(io,file=fname,form='formatted',status=stat_int,action='write',iostat=error,access=access_int)

if (error.ne.0) write(output_unit,'(3a,/)') 'The file ', fname,' could not be opened for writting'
rewind(io)

open_file_write=io

end function open_file_write

! function that open a file and return the unit for reading
integer function open_file_read(fname)
implicit none
character(len=*), intent(in) :: fname
! internal variables
logical :: i_file
integer :: io,error

open_file_read=-10
io=-1
i_file=inquire_file(fname)

if (i_file) then
! find unit of the file that is alread
  io=inquire_file_open(fname)
else
  return
endif
! find a new unit for the file

if (io.lt.0) io=inquire_file_unit(io)


open (io,file=fname,form='formatted',status='old',action='read',iostat=error)

if (error.ne.0) write(output_unit,'(3a,/)') 'The file ', fname,' could ne be opened for reading'
rewind(io)

open_file_read=io

end function open_file_read

! function that closes a file
subroutine close_file(filename,io)
implicit none
integer, intent(in) :: io
character(len=*), intent(in) :: filename
! internal
integer :: error

close(io,iostat=error)

if (error.ne.0) write(output_unit,'(/,3a,/)') 'The file ', filename, ' could not be closed'
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!
! PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!

! Inquire the presence of absence of a file in the working directory
logical function inquire_file(filename)
implicit none
character(len=*), intent(in) :: filename
! internal variables

inquire_file=.false.

inquire(file=filename,exist=inquire_file)

if (.not.inquire_file) then
  write(output_unit,'(/,3a)') 'file ', filename, ' has not been found'
endif

end function inquire_file

! Inquire if a file is open or not
integer function inquire_file_open(filename)
implicit none
character(len=*), intent(in) :: filename
! internal variables
logical :: i_open

i_open=.false.
inquire_file_open=-1

inquire(file=filename,opened=i_open,number=inquire_file_open)

end function inquire_file_open

! Inquire if the unit of a file and find a unit number if necessary
integer function inquire_file_unit(io)
implicit none
integer, intent(in) :: io
! internal variables
integer :: i,io_test
logical :: i_open

i_open=.True.
io_test=io


if (io_test.gt.0) inquire(unit=io_test,opened=i_open)

! if already open, find a new unit
if (i_open) then
  io_test=io+10
   do i=1,100
      io_test=io_test+1
      inquire(unit=io_test,opened=i_open)
      if (.not.i_open) exit
   enddo
endif

inquire_file_unit=io_test

end function inquire_file_unit

end module m_io_files_utils
