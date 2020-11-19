module m_topo_sd
use m_derived_types
use m_topo_commons
use m_convert
use m_io_utils
use m_io_files_utils

private
public :: get_charge_map
contains

subroutine get_charge_map(signature)
implicit none
integer, intent(in) :: signature
! internal
real(kind=8), allocatable, dimension(:,:) :: topo_map
integer :: shape_Q_topo(2),N_cell
integer :: i,io
character(len=50) :: name

shape_Q_topo=shape(Q_topo)
N_cell=shape_Q_topo(2)
allocate(topo_map(5,N_cell))
topo_map=0.0d0

do i=1,N_cell
   topo_map(:,i)=get_charge(i)
enddo

name=convert('topomap_',signature,'.dat')
io=open_file_write(name)

call dump_config(io,topo_map)

call close_file(name,io)

end subroutine get_charge_map

end module m_topo_sd
