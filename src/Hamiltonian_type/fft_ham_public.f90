module m_fft_H_public
use m_fft_H_base,  only: fft_H
use m_fft_H_fftw,  only: fft_H_fftw
#ifdef CPP_CUDA
use m_fft_H_cufft, only: fft_H_cufft
#endif
implicit none
public

integer,private ::  mode=-1 !stores which implementation is used (1:FFTW3, 2: cuFFT, 0: none, -1: not initialized)

interface  get_fft_H
    module procedure get_fft_H_single
    module procedure get_fft_H_N
end interface
private :: get_fft_H_single, get_fft_H_N

contains

    subroutine set_fft_H_mode_io(fname_in)
        use m_io_files_utils, only : open_file_read,close_file
        use m_io_utils,only: get_parameter
        character(*),intent(in),optional    :: fname_in
        !internal
        character(*),parameter              :: fname_default='input'
        character(:), allocatable           :: fname

        integer ::  io_param

#ifdef CPP_CUDA
        mode=2
#elif defined CPP_FFTW3
        mode=1
#else
        mode=0
#endif
        if(present(fname_in))then
            fname=fname_in
        else
            fname=fname_default
        endif
        io_param=open_file_read(fname)
        call get_parameter(io_param,fname,'Hamiltonian_fft_mode',mode)

        call close_file(fname,io_param)

        select case(mode)
        case(1)
            write(*,'(/A/)') "Using Hamiltonian_fft implementation: FFTW3"
        case(2)
            write(*,'(/A/)') "Using Hamiltonian_fft implementation: Cuda"
        case(0)
            write(*,'(/A/)') "No implementation for Hamiltonian_fft mode available"
            !no implementation is fine as long as it is not actually allocated
            continue
        case default
            ERROR STOP "READ IN UNIMPLEMENTED Hamiltonian_fft_mode"
        end select
    end subroutine

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
        case(0)
            ERROR STOP "Cannot allocate Fourier Hamiltonian type, requires either CPP_FFTW3 or CPP_CUDA"
        case(-1)
            ERROR STOP "Cannot allocate Fourier Hamiltonian type, mode has not been set using set_fft_H_mode_io"
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
        case(0)
            ERROR STOP "Cannot allocate Fourier Hamiltonian type, requires either CPP_FFTW3 or CPP_CUDA"
        case(-1)
            ERROR STOP "Cannot allocate Fourier Hamiltonian type, mode has not been set using set_fft_H_mode_io"
        case default
            ERROR STOP "UNEXPECTED MODE SET"
        end select

    end subroutine

    subroutine bcast_fft_H_mode(com)
        use mpi_util
        class(mpi_type),intent(in)  :: com

        Call bcast(mode,com)
    end subroutine


end module
