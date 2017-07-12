      subroutine MC_transfer(Ilat,Spin_store,irank_discuss,spin,shape_spin)
      use m_mpi_prop, only : MPI_COMM_BOX
      implicit none
! position of the spin
      integer, intent(in) :: Ilat(4),shape_spin(5)
! value of the spin that will be send
      real(kind=8), intent(in) :: Spin_store(3)
! the spin is sent from
      integer, intent(in) :: irank_discuss
! spin which has to be place into the mnatrix of spins
       real(kind=8), intent(inout) :: spin(shape_spin(1),shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5))
! internal variables
       integer :: ierr,Ilat_local(4)
       real(kind=8) :: S_trans(3)


      include 'mpif.h'

      Ilat_local=Ilat
      S_trans=Spin_store
! the rank 0 of irank_nei must be the domain in the center and the neighbors have the rank 1,2,3...
! all communications are none blocking.
      call mpi_bcast(Ilat_local,4,MPI_INTEGER,irank_discuss,MPI_COMM_BOX,ierr)
      call mpi_bcast(S_trans,3,MPI_REAL8,irank_discuss,MPI_COMM_BOX,ierr)

! update the spin
      spin(4:6,Ilat_local(1),Ilat_local(2),Ilat_local(3),Ilat_local(4))=S_trans

       end subroutine MC_transfer
