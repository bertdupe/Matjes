module m_mpi_start_end
use mpi_basic

contains

subroutine init_MPI(mpi_world)
    type(mpi_type),intent(inout)  :: mpi_world

#ifdef CPP_MPI
    integer         :: provided
    integer         :: ierr

!    Call MPI_init(ierr)
    Call MPI_INIT_THREAD(MPI_THREAD_FUNNELED, provided, ierr)
    if(ierr/=0) STOP "FAILED initialize MPI"
    if(provided/=MPI_THREAD_FUNNELED) STOP "failed to initialize MPI with MPI_THREAD_FUNNELED"
    Call mpi_world%set(MPI_COMM_WORLD)
#endif
end subroutine


subroutine end_MPI()
    integer         :: ierr

#ifdef CPP_MPI
    Call MPI_finalize(ierr)
    if(ierr/=0) STOP "Failed finalize mpi"
#else
    continue
#endif


end subroutine

end module
