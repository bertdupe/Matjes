module m_FFTW3
    use, intrinsic :: iso_c_binding
#ifdef CPP_FFTW3
    include 'fftw3.f03'
#endif
    public

contains 
subroutine fftw_init()
#ifdef CPP_FFTW3
!$  integer :: ierr
!$  ierr=fftw_init_threads()
!$  if(ierr==0) STOP "ERROR WITH OPENMP INITIALIZATION OF FFTW3"
#else
    continue
#endif
end subroutine
end module
