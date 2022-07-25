module m_fftw3
    use, intrinsic :: iso_c_binding
#if defined(CPP_FFTW3) && !defined(CPP_FFTWMPI)
    include 'fftw3.f03'
#endif
#ifdef CPP_FFTWMPI
    include 'fftw3-mpi.f03'
#endif
    public

contains 

subroutine fftw_init()
    integer :: ierr

#if defined(CPP_FFTWMPI)
    call fftw_mpi_init()
#endif

#if defined(CPP_FFTW3_THREAD)
!$  ierr=fftw_init_threads()
!$  if(ierr==0) STOP "ERROR WITH OPENMP INITIALIZATION OF FFTW3"
#endif

    continue

end subroutine

end module
