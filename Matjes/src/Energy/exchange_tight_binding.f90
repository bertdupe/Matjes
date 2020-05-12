module m_exchange_TB
use m_symmetry_operators
use m_Hamiltonian_variables, only : coeff_ham_inter_spec

type(coeff_ham_inter_spec), target, public, protected :: TB

private
public :: get_ham_TB

contains

subroutine get_ham_TB(fname)
use m_io_utils
use m_io_files_utils
use m_convert
implicit none
character(len=*), intent(in) :: fname
! internal
integer :: io_param,neighbor_hoping
integer :: i
real(kind=8), allocatable :: t_local(:)

TB%c_ham=1.0d0
TB%name='Tight-binding'
TB%N_shell=-1
TB%order=2

io_param=open_file_read(fname)
call get_parameter(io_param,fname,'c_tij',TB%c_ham)
! count the exchange coefficients if present
! count the number of anisotropy coefficients if present
call get_parameter(io_param,fname,'N_hoping',TB%N_shell)
if (TB%N_shell.eq.-1) then
   neighbor_hoping=count_variables(io_param,'t_',fname)
else
   neighbor_hoping=TB%N_shell
endif


if (neighbor_hoping.ne.0) then
   allocate(t_local(neighbor_hoping))
   t_local=0.0d0

   call get_coeff(io_param,fname,'t_',t_local)
   neighbor_hoping=number_nonzero_coeff(t_local,'TB')
endif
if (neighbor_hoping.ne.0) TB%i_exist=.true.

!
! allocate the number of shell in TB
!
allocate(TB%ham(neighbor_hoping))
do i=1,neighbor_hoping
   allocate(TB%ham(i)%H(2,2))
   TB%ham(i)%H=0.0d0
enddo

! get the symmetric exchange
if (neighbor_hoping.ne.0) then
  do i=1,neighbor_hoping
     TB%ham(i)%H(1,2)=TB%c_ham*t_local(i)
     TB%ham(i)%H(2,1)=TB%c_ham*t_local(i)
  enddo
endif

end subroutine get_ham_TB

end module m_exchange_TB
