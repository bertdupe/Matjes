module m_initialize_path
!use m_gneb_parameters
use m_io_gneb
use m_type_lattice,only: lattice
use m_derived_types, only : io_parameter
use m_rotation, only : rotation_axis
use m_convert
use m_path
use m_operator_pointer_utils
use m_write_spin
use m_createspinfile
use m_energyfield
use m_minimize
use m_H_public
use m_rw_GNEB,only: GNEB_input  
implicit none
private
public :: path_initialization

contains
subroutine path_initialization(images,io_simu,io_gneb,Ham)
    implicit none
    type(io_parameter), intent(in)  :: io_simu
    type(GNEB_input),intent(in)     :: io_gneb
    type(lattice), intent(inout)    :: images(:)
    class(t_H),intent(in)           :: Ham(:)
    ! internal variable
    integer :: i_nim,nim,N_cell
    logical :: exists=.False.,found=.False.
    character(35) :: num
    character(:),allocatable        ::  momfile_i,momfile_f
    
    nim=size(images)
    select case(io_gneb%initpath)
     case(3)
!PB: WHY DOES THIS FIRST TRY TO READ IN ALL FILES?
       write (6,'(a)') "Initial guess for the path from file"
       call read_path(io_gneb%restartfile_path,images,exists)
       if (.not.exists) stop
       momfile_i = convert(io_gneb%restartfile_path,'_1.dat')
       momfile_f = convert(io_gneb%restartfile_path,'_',nim,'.dat')
       call read_inifin(momfile_i,momfile_f,images)
#if 0    
       call WriteSpinAndCorrFile(path(:,:,1),'SpinSTM_GNEB_ini.dat')
       call CreateSpinFile(path(:,:,1),'povray_GNEB_ini.dat')
       call WriteSpinAndCorrFile(path(:,:,nim),'SpinSTM_GNEB_fin.dat')
       call CreateSpinFile(path(:,:,nim),'povray_GNEB_fin.dat')
    
       write (6,'(a)') "Relaxing the first image via the infinite damping method..."
    
       call minimize_infdamp(path(:,:,1),io_simu,Ham)
       write (6,'(a)') "Done!"
    
       write (6,'(a)') "Relaxing the first image via the infinite damping method..."
       call minimize_infdamp(path(:,:,nim),io_simu,Ham)
       write (6,'(a)') "Done!"
    
     case(2)
    
       write (6,'(a)') "Initial guess for the path from file"
       call read_path(restartfile_path,path,exists)
       if (.not.exists) stop
       momfile_i = convert(restartfile_path,'_1.dat')
       momfile_f = convert(restartfile_path,'_',nim,'.dat')
       call read_inifin(momfile_i,momfile_f,path)
       call WriteSpinAndCorrFile(path(:,:,1),'SpinSTM_GNEB_ini.dat')
       call CreateSpinFile(path(:,:,1),'povray_GNEB_ini.dat')
       call WriteSpinAndCorrFile(path(:,:,nim),'SpinSTM_GNEB_fin.dat')
       call CreateSpinFile(path(:,:,nim),'povray_GNEB_fin.dat')
    
       write (6,'(a)') "Relaxing the first image..."
    
       call minimize(path(:,:,1),io_simu,Ham)
       write (6,'(a)') "Done!"
    
       write (6,'(a)') "Relaxing the last image..."
       call minimize(path(:,:,nim),io_simu,Ham)
       write (6,'(a)') "Done!"
    
    
     case(1)
    
       write (6,'(a)') "Initial and final states are read in file"
       write (6,'(a)') "no minimization of the energy"
    
       call read_inifin(momfile_i,momfile_f,path)
       call WriteSpinAndCorrFile(path(:,:,1),'SpinSTM_GNEB_ini.dat')
       call CreateSpinFile(path(:,:,1),'povray_GNEB_ini.dat')
       call WriteSpinAndCorrFile(path(:,:,nim),'SpinSTM_GNEB_fin.dat')
       call CreateSpinFile(path(:,:,nim),'povray_GNEB_fin.dat')
    
       call geodesic_path(amp_rnd_path,path)
       write (*,'(a)') "Geodesic path generated"
    
     case default
    
       write (6,'(a)') "Initial and final states are read in file"
       write (6,'(a)') "Energy is minimized before GNEB"
    
       call read_inifin(momfile_i,momfile_f,path)
       call WriteSpinAndCorrFile(path(:,:,1),'SpinSTM_GNEB_ini.dat')
       call CreateSpinFile(path(:,:,1),'povray_GNEB_ini.dat')
       call WriteSpinAndCorrFile(path(:,:,nim),'SpinSTM_GNEB_fin.dat')
       call CreateSpinFile(path(:,:,nim),'povray_GNEB_fin.dat')
    
       write (*,'(a)') "Relaxing the first image..."
    
       call minimize(path(:,:,1),io_simu,Ham)
       write (*,'(a)') "Done!"
    
       write (*,'(a)') "Relaxing the last image..."
       call minimize(path(:,:,nim),io_simu,Ham)
       write (*,'(a)') "Done!"
    
       call geodesic_path(amp_rnd_path,path)
       write (*,'(a)') "Geodesic path generated"
#endif
    end select

    ERROR STOP "FINISH THIS"

end subroutine path_initialization

end module m_initialize_path
