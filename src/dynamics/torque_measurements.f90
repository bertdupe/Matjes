module m_torque_measurements
use m_io_files_utils
use m_io_utils

integer :: io_convergence
contains

subroutine init_Torque_measure()

   ! file for convergence criteria
   io_convergence=open_file_write('convergence.dat')

end subroutine



subroutine get_Torque_measure(j,Edy,Torque)
real(8), intent(in)                 :: Edy
integer, intent(in)                 :: j
real(8),dimension(:,:),intent(in)   :: Torque

integer :: N
real(8) :: ave_torque,max_torque

   N=size(Torque,2)
   ave_torque=(sum(Torque**2))/real(N,8)
   max_torque=maxval(Torque**2)
   write(io_convergence,'(I10,3x,3(E20.12E3,3x))') j,Edy,max_torque,ave_torque

end subroutine





subroutine finalize_Torque_measure()

   call close_file('convergence.dat',io_convergence)

end subroutine


end module m_torque_measurements
