subroutine rw_lattice(my_lattice)
use m_vector, only : norm,cross
use m_constants, only : pi
use m_derived_types
use m_io_files_utils
use m_io_utils
implicit none
! out
type(lattice), intent(out) :: my_lattice
! dummy variable
integer :: j,dimension(3),io
integer :: N_site
logical :: Periodic_log(3),fexist
real (kind=8) :: a(3),ss,net(3,3)

write(6,'(a)') 'reading the lattice parameters in the input file'
inquire(file='input',exist=fexist)
if(.not. fexist)then
    write(*,*) "Did not find input file, you should at least supply some lattices"
    STOP "Supply input file"
endif
io=open_file_read('input')

call get_parameter(io,'input','Periodic_log',3,Periodic_log)
call get_parameter(io,'input','alat',3,a)

call get_parameter(io,'input','lattice',3,3,net)
do j=1,3
   ss=norm(net(j,:))
   if (ss.lt.1.0d-8) stop 'lattice vectors are crap'
   net(j,:)=net(j,:)/ss
enddo

call get_parameter(io,'input','Nsize',3,dimension)
if (dimension(1)*dimension(2)*dimension(3).eq.0) then
   write(6,*) 'dimension of lattice must be >0 along all directions'
   stop
elseif((dimension(1).ne.dimension(2)).and.(dimension(2).ne.1)) then
   write(6,*) 'asymetric cell choosen'
elseif((dimension(1).ne.dimension(2)).and.(dimension(2).eq.1)) then
   write(6,*) 'chain geometry choosen'
endif
N_site=dimension(1)*dimension(2)*dimension(3)

call close_file('input',io)


!! put the magnetic lattice in order
Call my_lattice%init_geo(net,a,dimension,Periodic_log)

end subroutine rw_lattice
