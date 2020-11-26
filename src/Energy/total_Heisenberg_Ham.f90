module m_total_Heisenberg_Ham
#if 0
use m_summer_exp
use m_Hamiltonian_variables, only : coeff_ham_inter_spec

type(coeff_ham_inter_spec), target, public :: ham_tot_heisenberg

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
#endif
end module m_total_Heisenberg_Ham
