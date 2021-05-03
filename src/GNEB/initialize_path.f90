module m_initialize_path
!use m_gneb_parameters
use m_io_gneb
use m_type_lattice,only: lattice
use m_derived_types, only : io_parameter
use m_rotation, only : rotation_axis
use m_convert
use m_path
use m_write_spin
use m_createspinfile
use m_energyfield
use m_minimize
use m_H_public
use m_rw_GNEB,only: GNEB_input  
use m_hamiltonian_collection, only: hamiltonian
implicit none
private
public :: path_initialization

contains
subroutine path_initialization(images,io_simu,io_gneb,H)
    implicit none
    type(io_parameter), intent(in)  :: io_simu
    type(GNEB_input),intent(in)     :: io_gneb
    type(lattice), intent(inout)    :: images(:)
    type(hamiltonian),intent(inout) :: H
    ! internal variable
    integer :: nim
    
    nim=size(images)

    if(io_gneb%read_path)then
        call read_path(io_gneb%restartfile_path,images)
    endif
    if(io_gneb%read_outer)then
        Call read_path_inifin(io_gneb,images)
    endif

    call WriteSpinAndCorrFile(images(1)%M%modes_v,'SpinSTM_GNEB_ini.dat')
    call CreateSpinFile(images(1)%M%modes_v,'povray_GNEB_ini.dat')
    call WriteSpinAndCorrFile(images(nim)%M%modes_v,'SpinSTM_GNEB_fin.dat')
    call CreateSpinFile(images(nim)%M%modes_v,'povray_GNEB_fin.dat')

    select case(io_gneb%min_type)
    case(0)
       write(*,*) "WARNING, NO INITIAL MINIMIZATION FOR GNEB CHOSEN"
    case(1)
        write (6,'(a)') "Relaxing the first image via the infinite damping method..."
        call minimize_infdamp_run(images(1),io_simu,io_gneb%io_min,H)
        write (6,'(a)') "Done!"
        write (6,'(a)') "Relaxing the first image via the infinite damping method..."
        call minimize_infdamp_run(images(nim),io_simu,io_gneb%io_min,H)
        write (6,'(a)') "Done!"
    case(2)
        write (6,'(a)') "Relaxing the first image..."
        call minimize_run(images(1),io_simu,io_gneb%io_min,H)
        write (6,'(a)') "Done!"
        write (6,'(a)') "Relaxing the last image..."
        call minimize_run(images(nim),io_simu,io_gneb%io_min,H)
        write (6,'(a)') "Done!"
    case default
        ERROR STOP "UNEXPECTED io_gneb%min_type"
    end select

    if(.not.io_gneb%read_path)then
        call geodesic_path(io_gneb%amp_rnd_path,images)
        write (*,'(a)') "Geodesic path generated"
    endif

end subroutine path_initialization

end module m_initialize_path
