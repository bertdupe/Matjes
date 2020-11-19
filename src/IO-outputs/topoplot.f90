module m_topoplot
use m_topo_commons
use m_convert
use m_io_files_utils

interface get_topoplot
  module procedure get_topoplot_int,get_topoplot_out
end interface

private
public :: get_topoplot

contains


subroutine get_topoplot_int(file,tag)
implicit none
integer, intent(in) :: tag
character(len=*), intent(in) :: file
! inteernal variables
integer :: i
integer :: shape_Q(2),io_out
character(len=50) :: fname

shape_Q=shape(Q_topo)
fname=convert(file,tag,'.dat')

io_out=open_file_write(fname)

do i=1,shape_Q(1)

   Write(io_out,'(5(2x,E20.12E3))') get_charge(i)

enddo

call close_file(fname,io_out)

end subroutine get_topoplot_int


subroutine get_topoplot_out(matrix)
implicit none
real(kind=8), intent(inout) :: matrix(:,:)
! inteernal variables
integer :: i
integer :: shape_matrix(2)

shape_matrix=shape(matrix)

do i=1,shape_matrix(2)

   matrix(:,i)=get_charge(i)

enddo

end subroutine get_topoplot_out

end module m_topoplot
