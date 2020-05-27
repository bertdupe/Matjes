module m_print_Beff
use m_io_files_utils
use m_convert

private
public :: print_Beff

contains

subroutine print_Beff(tag,B)
implicit none
integer, intent(in) :: tag
real(kind=8) , intent(in) :: B(:,:)
! internal
integer :: shape_B(2),io_file,iomp
character(len=50) :: file_name,form

shape_B=shape(B)
form=convert('(',shape_B(1),'(E20.12E3,2x))')
file_name=convert('Beff_',tag,'.dat')
io_file=open_file_write(file_name)

do iomp=1,shape_B(2)
   write(io_file,form) B(:,iomp)
enddo

call close_file(file_name,io_file)

end subroutine print_Beff

end module
