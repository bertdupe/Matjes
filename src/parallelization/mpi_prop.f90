      module m_mpi_prop
! MPI communicator used in the calculation
       integer :: MPI_COMM
! MPI group
       integer :: all_world,working_world,working_group
! rank and irank of the processes
       integer :: irank_working,isize_working,irank,isize
! setup for the cartesian communicator
       integer :: MPI_COMM_CART
       integer, allocatable :: coords(:)
! size of the domains
       integer :: width,length,height
! mpi group and comm for each box
      integer :: working_box, MPI_COMM_BOX,irank_box
! mpi group and comm for neighbors
!      integer :: MPI_COMM_NEI,group_nei
! group that keep tracks of the master ranks
      integer :: MPI_COMM_MASTER,group_master
! size of the domain and starting point for the MC
      integer :: N(3), start(3)
      end module m_mpi_prop
