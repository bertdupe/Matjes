module m_Hamiltonian_variables
use m_basic_types

! Hamiltonian coefficients
! to be deleted
type Coeff_Ham
     real(kind=8) :: c_Ji=-1.0d0
     real(kind=8) :: c_DM=-1.0d0
     real(kind=8) :: c_JB=-1.0d0
     real(kind=8) :: c_Ki=-1.0d0
     real(kind=8) :: c_ani=1.0d0
! Exchange interaction
     real(kind=8), allocatable, dimension(:,:,:) :: exchange
! DMI interaction
     real(kind=8), allocatable, dimension(:,:,:) :: DMI
! magnetocrystalline interaction
     real(kind=8), allocatable, dimension(:,:,:) :: ani
! magnetocrystalline interaction
     real(kind=8), allocatable, dimension(:,:) :: Zeeman
! total Hamiltonian
     type(shell_Ham), allocatable, dimension(:) :: total_shell
! stoner parameter
     real(kind=8) :: Ist=0.0d0
! biquadratic interaction
     real(kind=8) :: Biq=0.0d0
! 4-spin interaction
     real(kind=8) :: fours=0.0d0
! presence or absence of interactions
     logical :: i_DM=.false.
     logical :: i_four=.false.
     logical :: i_biq=.false.
     logical :: i_dip=.false.
     logical :: i_exch=.false.
     logical :: i_ani=.false.
end type Coeff_Ham

!!!!!!!!
! Hamiltonian type
!!!!!!!!

type coeff_ham_inter_spec
     real(kind=8) :: c_ham=-1.0d0
     integer :: N_shell,order
     character(len=30) :: name=''
     logical :: i_exist=.false.
     type(site_Ham), allocatable, dimension(:) :: ham
end type coeff_ham_inter_spec

type coeff_ham_inter_spec_pointer
     character(len=30) :: name=''
     integer :: order
     logical :: onsite=.false.
     type(Op_real), allocatable :: ham(:)
end type coeff_ham_inter_spec_pointer


!
! number of Hamiltonian in one shell
!
type shell_Ham
     type(site_Ham), dimension(:), allocatable :: atom
     ! size of the Hamiltonian of order N
     integer :: line,column
end type

type shell_Ham_order_N
     type(shell_Ham), dimension(:), allocatable :: order
     ! number of order in the Hamiltonian
     integer, allocatable,dimension(:) :: num
end type shell_Ham_order_N

type H_vois
     type(shell_Ham_order_N), dimension(:), allocatable :: shell_num
     ! number of bounds/shell
     integer, dimension(:), allocatable :: num
end type

end module m_Hamiltonian_variables
