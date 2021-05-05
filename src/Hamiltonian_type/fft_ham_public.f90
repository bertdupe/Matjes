module m_fft_H_public
use m_fft_H_base, only: fft_H
use m_fft_H_fftw, only: fft_H_fftw
implicit none
public

interface  get_fft_H
    module procedure get_fft_H_single
    module procedure get_fft_H_N
end interface
private :: get_fft_H_single, get_fft_H_N


contains
    subroutine get_fft_H_single(H_out)
        class(fft_H),intent(out),allocatable      :: H_out 

        allocate(fft_H_fftw::H_out)
    end subroutine
    
    subroutine get_fft_H_N(H_out,N)
        class(fft_H),intent(out),allocatable    :: H_out (:)
        integer,intent(in)                      :: N
        allocate(fft_H_fftw::H_out(N))
    end subroutine

end module
