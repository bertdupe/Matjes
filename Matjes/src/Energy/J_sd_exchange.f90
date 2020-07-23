module m_J_sd_exchange

use m_symmetry_operators
use m_Hamiltonian_variables, only : coeff_ham_inter_spec
use m_rw_TB
use m_lattice, only : my_order_parameters

type(coeff_ham_inter_spec), target, public, protected :: Jsd_ham_TB
complex(kind=8), allocatable, dimension(:,:) , public, protected :: Jsd_mat

private
public :: get_Jsd_ham_TB

contains

subroutine get_Jsd_ham_TB(fname,dim_ham)
use m_convert
implicit none
integer, intent(in) :: dim_ham
character(len=*), intent(in) :: fname

! internal
character(len=50) :: form
integer :: i
integer :: x_start, x_end, y_start, y_end

x_start=-1
x_end=-1
y_start=-1
y_end=-1

! Multiplicative coefficient
Jsd_ham_TB%c_ham=1.0d0
! Name of the model
Jsd_ham_TB%name='J_sd-Exchange'
Jsd_ham_TB%N_shell=-1
! Order of the Hamiltonian (tensor of rank "order")
Jsd_ham_TB%order=2

if (.not. allocated(TB_params%Jsd)) return
if (count(abs(TB_params%Jsd).gt.1.0d-8).eq.0) return

Jsd_ham_TB%i_exist=.true.

call get_borders('magnetic',x_start,x_end,'Tight-binding',y_start,y_end,my_order_parameters)

! Allocate the different blocs in the total Hamiltonian
allocate(Jsd_ham_TB%ham(1))
allocate(Jsd_ham_TB%ham(1)%H(dim_ham,dim_ham**2),Jsd_mat(dim_ham,dim_ham**2))
Jsd_ham_TB%ham(1)%H=0.0d0
Jsd_mat=cmplx(0.0d0)

! in this form the J_sd coupling is a rank 3 tensor

do i=1,size(TB_params%Jsd)
  Jsd_mat(x_start+2,3*dim_ham+y_start)=TB_params%Jsd(i)
  Jsd_mat(x_start,3*dim_ham+y_start+1)=TB_params%Jsd(i)
  Jsd_mat(x_start+1,3*dim_ham+y_start+1)=-cmplx(0.0d0,TB_params%Jsd(i))

  Jsd_mat(x_start,4*dim_ham+y_start)=TB_params%Jsd(i)
  Jsd_mat(x_start+1,4*dim_ham+y_start)=cmplx(0.0d0,TB_params%Jsd(i))
  Jsd_mat(x_start+2,4*dim_ham+y_start+1)=-TB_params%Jsd(i)
enddo

form=convert('(',dim_ham,'(f12.8,2x))')
write(6,'(a)') ''
write(6,'(a)') 'J_sd tight-binding (real part)'
do i=1,dim_ham**2
   write(6,form) real(Jsd_mat(:,i))
enddo
write(6,'(a)') 'J_sd tight-binding (complex part)'
do i=1,dim_ham**2
   write(6,form) aimag(Jsd_mat(:,i))
enddo
write(6,'(a)') ''

end subroutine get_Jsd_ham_TB

end module m_J_sd_exchange
