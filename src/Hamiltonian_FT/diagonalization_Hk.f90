module m_diagonalization_Hk
use m_parameters_rw_high
use m_set_Hamiltonian_FT
!
! 1. make the FT of each of the Ham individually
!  - create one type per Ham
!  - collect the real space differences for each of them
! 2. combine all the Hamiltonian and real space distances and the into a big Hamiltonian
! 3. for each of the k, get the phase factor and put it into the Hamiltonian
! 4. diagonalise Hk
!
use m_construction_Hk
private
public  :: diagonalize_Ham

contains

subroutine diagonalize_Ham()

   call set_Hamiltonians(Ham_res,Ham_comb,keep_res,H_io,lat)
   write(*,*) 'toto'
   stop

end subroutine

end module m_diagonalization_Hk
