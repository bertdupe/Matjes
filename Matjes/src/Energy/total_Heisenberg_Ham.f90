module m_total_Heisenberg_Ham
use m_exchange_heisenberg
use m_anisotropy_heisenberg
use m_zeeman
use m_anisotropy_heisenberg
use m_summer_exp
use m_Hamiltonian_variables, only : coeff_ham_inter_spec

type(coeff_ham_inter_spec), target, public, protected :: ham_tot_heisenberg

private
public :: get_total_Heisenberg_Ham

contains

subroutine get_total_Heisenberg_Ham(fname,dim_ham,Ms)
use m_convert
implicit none
character(len=*), intent(in) :: fname
integer, intent(in) :: dim_ham
real(kind=8), intent(in) :: Ms
! Internal variables
integer :: i,j
integer :: n_shell
character(len=50) :: form


ham_tot_heisenberg%order=2
ham_tot_heisenberg%name='heisenberg'

! get the anisotropy Hamiltonian
call get_ham_anisotropy(fname,dim_ham)
! get the zeeman
call get_ham_zeeman(fname,dim_ham,Ms)
! get the exchange Hamiltonian
call get_ham_exchange(fname,dim_ham)

n_shell=size(exchange%ham)+1

if (n_shell.ne.0) then
  ham_tot_heisenberg%i_exist=.true.
  ham_tot_heisenberg%N_shell=n_shell
else
  write(6,'(/a/)') 'no Heisenberg Hamiltonian found'
  return
endif

!!!!
!!!! The Heisenberg Hamiltonian exists
!!!!

allocate(ham_tot_heisenberg%ham(n_shell))
do i=1,n_shell
   allocate(ham_tot_heisenberg%ham(i)%H(dim_ham,dim_ham))
   ham_tot_heisenberg%ham(i)%H=0.0d0
enddo

if (anisotropy%i_exist) ham_tot_heisenberg%ham(1)%H=ham_tot_heisenberg%ham(1)%H+anisotropy%ham(1)%H
if (Zeeman%i_exist) ham_tot_heisenberg%ham(1)%H=ham_tot_heisenberg%ham(1)%H+zeeman%ham(1)%H

if (exchange%i_exist) then
  do i=2,n_shell
    ham_tot_heisenberg%ham(i)%H=ham_tot_heisenberg%ham(i)%H+exchange%ham(i-1)%H
  enddo
endif

if (temperature_strasbourg%i_exist) then
  do i=1,n_shell
    ham_tot_heisenberg%ham(i)%H=ham_tot_heisenberg%ham(i)%H+temperature_strasbourg%ham(i)%H
  enddo
endif

form=convert('(',dim_ham,'(f12.8,2x))')
write(6,'(a)') ''
write(6,'(a)') 'Total Hamiltonian Heisenberg'
do j=1,size(ham_tot_heisenberg%ham)
  write(6,'(a,I8)') 'Now in shell ', j-1
  do i=1,dim_ham
    write(6,form) ham_tot_heisenberg%ham(j)%H(:,i)
  enddo
enddo
write(6,'(a)') ''

end subroutine get_total_Heisenberg_Ham

end module m_total_Heisenberg_Ham
