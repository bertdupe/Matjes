module m_FFTW3
    use, intrinsic :: iso_c_binding
#ifdef CPP_FFTW3
    include 'fftw3.f03'
#endif
end module
