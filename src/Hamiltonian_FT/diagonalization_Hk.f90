module m_diagonalization_Hk
use m_parameters_rw_high
use m_set_Hamiltonian_FT
use m_FT_Ham_public
use m_FT_Ham_base
use m_io_files_utils
use m_parameters_FT_Ham
! the following module is used for the TB. You will find the module in the directory tight-binding
use m_highsym, only : set_highs_path,mv_kpts
use m_input_H_types
use m_derived_types
use m_FT_Ham_public
use m_FT_Ham_coo_rtok_base
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
    type(io_h),intent(in)               :: H_io
    type(lattice), intent(in)           :: lat

    ! internal Hamiltonians
    type(H_inp_real_to_k),allocatable :: FT_Ham(:)
    class(FT_Ham_base),allocatable    :: Hk
    type(parameters_FT_HAM_IO)        :: io_H_diag         ! parameters for the diagonalization
    ! high symmetry lines
    type(parameters_IO_HIGHS)         :: high_lines

    ! dummy variables
    real(8)   :: k(3)
    integer   :: io_input
    integer   :: i,n_kpts,n_eigen
    real(8), allocatable :: kpts(:,:)
    complex(8),allocatable  :: eigenvalues(:)      ! array containing the eigenvalues
    complex(8),allocatable  :: eigenvectors(:,:)   ! array containing the eigenvectors
    character(len=100)       :: form

    ! initialization
    k=0.0d0
    io_input=open_file_read('input')
    call io_H_diag%read_file(io_input,'input')

    ! read the high symmetry lines
    call high_lines%read_file('q',io_input,'input')
    call close_file('input',io_input)
    call set_highs_path(lat,high_lines)
    call mv_kpts(kpts)

    ! prepare the Hamiltonian based on the coo matrices for the FT
    call set_Hamiltonians_FT(FT_Ham,H_io,lat)     ! choose with which algoritm you want to work

    Call set_H(Hk,io_H_diag)   ! choose the Hamiltonian with which you would like to work (sparse, dense...)

    call Hk%init(FT_Ham,io_H_diag)    ! initialize the Hamiltonian matrix
    call Hk%set_work(eigenvalues,eigenvectors)

    io_input=open_file_write('dispersion.dat')
    write(form,'( "(3E20.12E3,", I10, "(x,E20.12E3,x,E20.12E3))" )') size(eigenvalues)

    n_kpts=size(kpts,2)
    do i=1,n_kpts
       call Hk%set_k(FT_Ham,kpts(:,i))
       call Hk%calc_eval(3,eigenvalues,n_eigen)
       write(io_input,form) kpts(:,i),real(eigenvalues),aimag(eigenvalues)
    enddo

    call close_file('dispersion.dat',io_input)
    stop 'dispersion done'
end subroutine

end module m_diagonalization_Hk
