module m_topo_sd
use m_derived_types
use m_topo_commons
use m_convert
use m_io_utils
use m_io_files_utils
implicit none

private
public :: get_charge_map
contains

subroutine get_charge_map(signature,lat,Q_neigh)
    implicit none
    integer, intent(in)         :: signature
    type(lattice), intent(in)   :: lat
    integer,intent(in)          :: Q_neigh(:,:) 
    
    ! internal
    real(kind=8), allocatable, dimension(:,:) :: topo_map
    integer :: N_cell
    integer :: i,io
    character(len=50) :: name
    
    N_cell=size(Q_neigh,2)
    allocate(topo_map(5,N_cell),source=0.0d0)
    
    do i=1,N_cell
       topo_map(:,i)=get_charge(i,lat,Q_neigh)
    enddo
    
    name=convert('topomap_',signature,'.dat')
    io=open_file_write(name)
    call dump_config(io,topo_map)
    call close_file(name,io)
    
end subroutine get_charge_map

end module m_topo_sd
