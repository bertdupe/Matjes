module m_diagonalization_Hk
use m_parameters_rw_high
use m_set_Hamiltonian_FT
use m_construction_Hk
use m_FT_Ham_base
use m_FT_Ham_coo
use m_io_files_utils
! the following module is used for the TB. You will find the module in the directory tight-binding
use m_highsym, only : set_highs_path,mv_kpts
implicit none

!
! 1. make the FT of each of the Ham individually
!  - create one type per Ham
!  - collect the real space differences for each of them
! 2. combine all the Hamiltonian and real space distances and the into a big Hamiltonian
! 3. for each of the k, get the phase factor and put it into the Hamiltonian
! 4. diagonalise Hk
!

private
public  :: diagonalize_Ham_FT

contains

subroutine diagonalize_Ham_FT(H_io,lat)
    use m_input_H_types
    use m_derived_types
    type(io_h),intent(in)               :: H_io
    type(lattice), intent(in)           :: lat

    ! internal Hamiltonians
    type(H_inp_real_to_k),allocatable :: FT_Ham(:)
!    type(H_inp_k_coo) :: FT_Ham_k

    ! high symmetry lines
    type(parameters_IO_HIGHS) :: high_lines

    ! dummy variables
    real(8)   :: k(3)
    integer   :: io_input
    integer   :: i
    real(8), allocatable :: kpts(:,:)

    ! initialization
    k=0.0d0

    ! read the high symmetry lines
    io_input=open_file_read('input')
    call high_lines%read_file('q',io_input,'input')
    call close_file('input',io_input)
    call set_highs_path(lat,high_lines)
    call mv_kpts(kpts)

    ! prepare the Hamiltonian for the FT
    call set_Hamiltonians_FT(FT_Ham,H_io,lat)

    ! get the phase of the Hamiltonian
    do i=1,size(kpts,2)
       call get_Hk(FT_Ham,kpts(:,i))
    enddo


end subroutine

end module m_diagonalization_Hk
