       module m_error
       interface error
#ifdef CPP_MPI
        module procedure err_comm
#endif
        module procedure err_alloc
       end interface error

       contains

#ifdef CPP_MPI
       subroutine err_comm(ierr,mpi_comm,COMM,routine)
       implicit none
       integer, intent(in) :: ierr,mpi_comm
       character, intent(in) :: routine,COMM
! internal variable

       write(6,'(a)') 'error in MPI comm'
       write(6,'(2a)') 'routine', COMM
       write(6,'(2a)') 'routine', routine
       write(6,'(a,I6)') 'error message', ierr

       end subroutine err_comm
#endif

       subroutine err_alloc(ierr, routine)
       implicit none
       integer, intent(in) :: ierr
       logical, intent(in) :: routine

       write(6,'(a)') 'error in allocation of memory'
       write(6,'(a)') 'routine', routine
       write(6,'(a,I6)') 'error message', ierr
       end subroutine err_alloc

       end module m_error
