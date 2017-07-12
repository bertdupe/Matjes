       module m_fft
       interface fft_exe
        module procedure fft_exe_r2c_1d
        module procedure fft_exe_r2c_3d
        module procedure fft_exe_c2r_3d
       end interface fft_exe
       contains

!!!
! part for the execution
!!!
       subroutine fft_exe_r2c_1d(N,rmat,cmat)
       implicit none
       real(kind=8), intent(in) :: rmat(:)
       complex, intent(inout) :: cmat(:)
       integer, intent(in) :: N
       !internal
       integer :: plan

       include 'fftw3.f'

       call dfftw_plan_dft_r2c_1d(plan,N,rmat,cmat,FFTW_ESTIMATE+FFTW_FORWARD)
       call dfftw_execute_dft_r2c(plan,rmat,cmat)

       end subroutine fft_exe_r2c_1d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

       function fft_exe_r2c_3d(plan,Nx,Ny,Nz,rmat)
       implicit none
       integer :: Nx,Ny,Nz
       real(kind=8) :: rmat(Nx,Ny,Nz)
       integer :: plan
       complex :: fft_exe_r2c_3d(Nx,Ny,Nz)

       include 'fftw3.f'

       call dfftw_execute_dft(plan,rmat,fft_exe_r2c_3d)

       end function fft_exe_r2c_3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       subroutine fft_exe_c2r_3d(rmat,Nx,Ny,Nz,cmat)
       implicit none
       integer, intent(in) :: Nx,Ny,Nz
       real(kind=8), intent(in) :: rmat(:,:,:)
       complex, intent(out) :: cmat(:,:,:)
! internal
       integer :: plan

       include 'fftw3.f'

       call dfftw_plan_dft_c2r_3d(plan,Nx,Ny,Nz,cmat,rmat,FFTW_ESTIMATE+FFTW_BACKWARD)
       call dfftw_execute_dft_c2r(plan,cmat,rmat)

       end subroutine fft_exe_c2r_3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       subroutine fft_exe_c2r_1d(rmat,Nx,cmat)
       implicit none
       integer, intent(in) :: Nx
       real(kind=8), intent(in) :: rmat(:)
       complex, intent(out) :: cmat(:)
! internal
       integer :: plan

       include 'fftw3.f'

       call dfftw_plan_dft_c2r_1d(plan,Nx,cmat,rmat,FFTW_ESTIMATE+FFTW_BACKWARD)
       call dfftw_execute_dft_c2r(plan,cmat,rmat)

       end subroutine fft_exe_c2r_1d

       end module m_fft
