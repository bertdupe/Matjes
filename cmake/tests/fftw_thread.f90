program test_fftw_thread
    USE, INTRINSIC :: ISO_C_BINDING 
    implicit none
    include 'fftw3.f03'

    integer     :: ierr

    ierr=fftw_init_threads()
    Call fftw_plan_with_nthreads(1)
end program
