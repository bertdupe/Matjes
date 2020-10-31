module m_get_table_nn
use m_derived_types, only: lattice
use m_io_files_utils, only: open_file_write,close_file
use m_io_utils, only: dump_config
use m_get_position, only : get_position
use m_mapping, only: mapping
use m_indexation, only: get_num_neighbors
use m_table_dist, only: get_table_of_distance
implicit none
contains
    subroutine get_table_nn(lat,N_Nneigh,indexNN_out,tableNN)
    ! get table_nn, so far only for tight-binding H creation
    ! done really fast and ugly
        integer,intent(in)                  :: N_Nneigh
        type(lattice),intent(in)            :: lat
        integer, allocatable,intent(out)    :: indexNN_out(:),tableNN(:,:,:,:,:,:)
        ! variable of the system
        real(kind=8), allocatable :: tabledist(:,:)
        real (kind=8), allocatable :: pos(:,:,:,:,:)
        integer ::  alloc_check
        integer :: tot_N_Nneigh
        integer, allocatable    :: indexNN(:,:)

        ! get the table of neighbors and the numbers of atom per shell
        allocate(tabledist(N_Nneigh,1),stat=alloc_check)
        if (alloc_check.ne.0) write(6,'(a)') 'out of memory cannot allocate table of distance'
        tabledist=0.0d0
        allocate(indexNN(N_Nneigh,1),stat=alloc_check)
        if (alloc_check.ne.0) write(6,'(a)') 'out of memory cannot allocate indexNN'
        indexNN=0
        call get_table_of_distance(lat%areal,N_Nneigh,lat%world,lat%cell,lat%nmag,tabledist)
        call get_num_neighbors(N_Nneigh,tabledist,lat%areal,lat%world,lat%cell,indexNN)
        
        ! allocate table of neighbors and masque
        tot_N_Nneigh=sum(indexNN(1:N_Nneigh,1),1)
        
        allocate(tableNN(5,tot_N_Nneigh,lat%dim_lat(1),lat%dim_lat(2),lat%dim_lat(3),lat%nmag),stat=alloc_check)
        if (alloc_check.ne.0) write(6,'(a)') 'out of memory cannot allocate tableNN'
        
        tableNN=0
        
        allocate(pos(3,lat%dim_lat(1),lat%dim_lat(2),lat%dim_lat(3),lat%nmag),stat=alloc_check)
        if (alloc_check.ne.0) write(6,'(a)') 'out of memory cannot allocate position for the mapping procedure'
        
        ! get position
        pos=0.0d0
        call get_position(pos,lat%dim_lat,lat%areal,lat%cell)
        call mapping(tabledist,N_Nneigh,lat%cell,indexNN,tableNN,lat)
        allocate(indexNN_out,source=indexNN(:,1))
    end subroutine
end module
