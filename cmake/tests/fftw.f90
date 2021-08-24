
program test_fftw
    USE, INTRINSIC :: ISO_C_BINDING 
    implicit none
    include 'fftw3.f03'

    integer(C_int)  :: howmany
    real(C_DOUBLE)              :: arr_in(3)
    complex(C_DOUBLE_COMPLEX)   :: arr_out(3)
    type(c_ptr)     :: plan

    arr_in=1.0d0
    arr_out=(2.0d0,1.0d0)
    howmany=int(3,C_INT)
    plan= fftw_plan_dft_r2c_1d(howmany, arr_in, arr_out,FFTW_FORWARD+FFTW_ESTIMATE)
    Call fftw_destroy_plan(plan)
end program
