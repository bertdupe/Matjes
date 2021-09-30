module m_user_info
implicit none

interface user_info
       module procedure user_info_raw
end interface user_info

private
public :: welcome,user_info

contains
! simple module to write welcome messages for the different parts of the code

subroutine welcome_message()
    use , intrinsic :: iso_fortran_env, only: output_unit
implicit none

    write(6,'(/,a)') 'Bonjour!! you are now using the Matjes code'
    write(6,'(a)') 'All the developers hope that you will enjoy this moment.'
    write(6,'(a)') 'If you have a problem, if you are happy with the code or if you wish to chat a bit'
    write(6,'(a)') 'send an email to bertrand.dupe@uliege.be'
#ifdef CPP_VERSIONGIT
    write(6,'(/2a)') "You are using the git-version: ",CPP_VERSIONGIT
#endif
    write(6,'(/)') 
end subroutine

subroutine welcome
    use , intrinsic :: iso_fortran_env, only:  compiler_options
#if CPP_MPI
use mpi_basic,only: mpi_world
    if(mpi_world%ismas)then
        Call welcome_message()
        write(6,'(a,I10,2x,a)') 'Congratulations, you are now using',mpi_world%NP,'MPI threads'
        Call print_omp(6)
        Call print_preprocessor_flags(6)
        write(6,'(A/A/)') "Used compiler options:", compiler_options()
    endif
#else
    Call welcome_message()
    Call print_omp(6)
    Call print_preprocessor_flags(6)
    write(6,'(A/A/)') "Used compiler options:", compiler_options()
#endif

end subroutine 

subroutine user_info_raw(io,time1,string,space)
implicit none
integer, intent(in) :: io
character(len=*), intent(in) :: string
logical, intent(in) :: space
real(kind=8), intent(inout) :: time1
! internal
character(len=30) :: forme
real(kind=8) :: time2

if (space) then
   forme='(/a/)'
else
   forme='(a)'
endif

if (time1.lt.1.0d-8) then
   call cpu_time(time1)
else
   call cpu_time(time2)
   write(io,'(a,E10.5,a)') 'duration: ',time2-time1,' seconds'
   time1=0.0d0
endif

!REINTRODUCE ONLY WRITING WITH RESPECTIVE MASTER... MIGHT REQUIRE MODULE INTERFACE OVERLOAD PROCEDURE WITH MPI_TYPE
!#ifdef CPP_MPI
!if (irank.eq.0)
!#endif
   write(io,forme) string
!#ifdef CPP_MPI
!endif
!#endif

end subroutine user_info_raw

subroutine print_omp(io_unit)
!$  use omp_lib, only: omp_get_max_threads
    integer,intent(in)  :: io_unit
!$      if(.true.)then
!$          write(6,'(A,I6,A/)') "Running with openMP using maximally ", omp_get_max_threads(), " threads."
!$      else
            write(6,'(A,/)') "Compiled without openMP."
!$      endif 
end subroutine

subroutine print_preprocessor_flags(io_unit)
    integer,intent(in)  :: io_unit
    character(len=*),parameter  :: forma='(2X,A)'

    write(io_unit,'(/A)') "Compiled using the following preprocessor flags:"
#ifdef CPP_BLAS
    write(io_unit,forma) "CPP_BLAS"
#endif
#ifdef CPP_CUDA
    write(io_unit,forma) "CPP_CUDA"
#endif
#ifdef CPP_DEBUG
    write(io_unit,forma) "CPP_DEBUG"
#endif
#ifdef CPP_EIGEN
    write(io_unit,forma) "CPP_EIGEN"
#endif
#ifdef CPP_FFTW3
    write(io_unit,forma) "CPP_FFTW3"
#endif
#ifdef CPP_FFTW3_THREAD
    write(io_unit,forma) "CPP_FFTW3_THREAD"
#endif
#ifdef CPP_MKL
    write(io_unit,forma) "CPP_MKL"
#endif
#ifdef CPP_MPI
    write(io_unit,forma) "CPP_MPI"
#endif
#ifdef CPP_MRG
    write(io_unit,forma) "CPP_MRG"
#endif
#ifdef CPP_SCRIPT
    write(io_unit,forma) "CPP_SCRIPT"
#endif
#ifdef CPP_USE_WORK
    write(io_unit,forma) "CPP_USE_WORK"
#endif
    write(io_unit,'(/)')

end subroutine


end module m_user_info
