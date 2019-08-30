      module m_randperm
      contains

      integer function randperm(irank,N,COM)
      implicit none
      integer, intent(in) :: N,irank,COM
! internal
      integer :: step,ierr(3)
      real(kind=8) :: Choice

      include 'mpif.h'

      CALL RANDOM_NUMBER(Choice)
      step=1+int(Choice*dble(N-1))
      call mpi_bcast(step,1,MPI_INTEGER,0,COM,ierr)

      call mpi_barrier(COM,ierr)

      randperm=mod(irank+step,N)

      end function randperm

      end module m_randperm
