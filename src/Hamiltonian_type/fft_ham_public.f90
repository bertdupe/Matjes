module m_fft_H_public
use m_fft_H_base,  only: fft_H
use m_fft_H_fftw,  only: fft_H_fftw
use m_fft_H_cufft, only: fft_H_cufft
implicit none
public

integer,private ::  mode=2  !chooses which implementation is used (1:FFTW3, 2: cuFFT)

interface  get_fft_H
    module procedure get_fft_H_single
    module procedure get_fft_H_N
end interface
private :: get_fft_H_single, get_fft_H_N


contains
    subroutine get_fft_H_single(H_out)
        class(fft_H),intent(out),allocatable      :: H_out 
        select case(mode)
        case(1)
#ifdef CPP_FFTW3
            allocate(fft_H_fftw::H_out)
#else
            ERROR STOP "CANNOT USE fft_H_fftw FOURIER IMPLEMENTATION WITHOUT FFTW (CPP_FFTW3)"
#endif
        case(2)
#ifdef CPP_CUDA
            allocate(fft_H_cufft::H_out)
#else
            ERROR STOP "CANNOT USE fft_H_fftw FOURIER IMPLEMENTATION WITHOUT FFTW (CPP_FFTW3)"
#endif
        case default
            ERROR STOP "UNEXPECTED MODE SET"
        end select
    end subroutine
    
    subroutine get_fft_H_N(H_out,N)
        class(fft_H),intent(out),allocatable    :: H_out (:)
        integer,intent(in)                      :: N

        select case(mode)
        case(1)
#ifdef CPP_FFTW3
            allocate(fft_H_fftw::H_out(N))
#else
            ERROR STOP "CANNOT USE fft_H_fftw FOURIER IMPLEMENTATION WITHOUT FFTW (CPP_FFTW3)"
#endif
        case(2)
#ifdef CPP_CUDA
            allocate(fft_H_cufft::H_out(N))
#else
            ERROR STOP "CANNOT USE fft_H_fftw FOURIER IMPLEMENTATION WITHOUT FFTW (CPP_FFTW3)"
#endif
        case default
            ERROR STOP "UNEXPECTED MODE SET"
        end select

    end subroutine

end module
