module m_diagonalization_Hk
use m_parameters_rw_high
use m_set_Hamiltonian_FT
use m_construction_Hk
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
    type(Hk_inp_t),allocatable :: FT_Ham(:)
    type(Hk_inp_t),allocatable :: FT_Ham_real(:),FT_Ham_complx(:)

    ! dummy variables
    real(8)                    :: k(3)

    ! initialization
    k=0.0d0


    ! prepare the Hamiltonian for the FT
    call set_Hamiltonians_FT(FT_Ham,H_io,lat)

    ! get the phase of the Hamiltonian
    call get_Hk(FT_Ham,k,FT_Ham_real,FT_Ham_complx)


end subroutine

end module m_diagonalization_Hk
