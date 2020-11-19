module m_initialize_path
use m_gneb_parameters
use m_io_gneb
use m_derived_types, only : io_parameter
use m_rotation, only : rotation_axis
use m_convert
use m_path
use m_operator_pointer_utils
use m_write_spin
use m_createspinfile
use m_energyfield
use m_minimize

contains

subroutine path_initialization(path,io_simu)
implicit none
type(io_parameter), intent(in) :: io_simu
real(kind=8), intent(inout) :: path(:,:,:)
! internal variable
integer :: shape_path(3)
integer :: i_nim,nim,N_cell
logical :: exists=.False.,found=.False.
character(35) :: num

shape_path=shape(path)
nim=shape_path(3)
N_cell=shape_path(2)

select case (initpath)

 case(3)

   write (6,'(a)') "Initial guess for the path from file"
   call read_path(restartfile_path,amp_rnd_path,path,exists)
   if (.not.exists) stop
   momfile_i = convert(restartfile_path,'_1.dat')
   momfile_f = convert(restartfile_path,'_',nim,'.dat')
   call read_inifin(momfile_i,momfile_f,amp_rnd,path)
   call WriteSpinAndCorrFile(path(:,:,1),'SpinSTM_GNEB_ini.dat')
   call CreateSpinFile(path(:,:,1),'povray_GNEB_ini.dat')
   call WriteSpinAndCorrFile(path(:,:,nim),'SpinSTM_GNEB_fin.dat')
   call CreateSpinFile(path(:,:,nim),'povray_GNEB_fin.dat')

   write (6,'(a)') "Relaxing the first image via the infinite damping method..."

   call minimize_infdamp(path(:,:,1),io_simu)
   write (6,'(a)') "Done!"

   write (6,'(a)') "Relaxing the first image via the infinite damping method..."
   call minimize_infdamp(path(:,:,nim),io_simu)
   write (6,'(a)') "Done!"


 case(2)

   write (6,'(a)') "Initial guess for the path from file"
   call read_path(restartfile_path,amp_rnd_path,path,exists)
   if (.not.exists) stop
   momfile_i = convert(restartfile_path,'_1.dat')
   momfile_f = convert(restartfile_path,'_',nim,'.dat')
   call read_inifin(momfile_i,momfile_f,amp_rnd,path)
   call WriteSpinAndCorrFile(path(:,:,1),'SpinSTM_GNEB_ini.dat')
   call CreateSpinFile(path(:,:,1),'povray_GNEB_ini.dat')
   call WriteSpinAndCorrFile(path(:,:,nim),'SpinSTM_GNEB_fin.dat')
   call CreateSpinFile(path(:,:,nim),'povray_GNEB_fin.dat')

   write (6,'(a)') "Relaxing the first image..."

   call minimize(path(:,:,1),io_simu)
   write (6,'(a)') "Done!"

   write (6,'(a)') "Relaxing the last image..."
   call minimize(path(:,:,nim),io_simu)
   write (6,'(a)') "Done!"


 case(1)

   write (6,'(a)') "Initial and final states are read in file"
   write (6,'(a)') "no minimization of the energy"

   call read_inifin(momfile_i,momfile_f,amp_rnd,path)
   call WriteSpinAndCorrFile(path(:,:,1),'SpinSTM_GNEB_ini.dat')
   call CreateSpinFile(path(:,:,1),'povray_GNEB_ini.dat')
   call WriteSpinAndCorrFile(path(:,:,nim),'SpinSTM_GNEB_fin.dat')
   call CreateSpinFile(path(:,:,nim),'povray_GNEB_fin.dat')

   call geodesic_path(amp_rnd_path,path)
   write (*,'(a)') "Geodesic path generated"

 case default

   write (6,'(a)') "Initial and final states are read in file"
   write (6,'(a)') "Energy is minimized before GNEB"

   call read_inifin(momfile_i,momfile_f,amp_rnd,path)
   call WriteSpinAndCorrFile(path(:,:,1),'SpinSTM_GNEB_ini.dat')
   call CreateSpinFile(path(:,:,1),'povray_GNEB_ini.dat')
   call WriteSpinAndCorrFile(path(:,:,nim),'SpinSTM_GNEB_fin.dat')
   call CreateSpinFile(path(:,:,nim),'povray_GNEB_fin.dat')

   write (*,'(a)') "Relaxing the first image..."

   call minimize(path(:,:,1),io_simu)
   write (*,'(a)') "Done!"

   write (*,'(a)') "Relaxing the last image..."
   call minimize(path(:,:,nim),io_simu)
   write (*,'(a)') "Done!"

   call geodesic_path(amp_rnd_path,path)
   write (*,'(a)') "Geodesic path generated"
end select

end subroutine path_initialization

end module m_initialize_path
