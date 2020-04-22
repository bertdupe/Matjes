module m_zeeman
use m_symmetry_operators
use m_lattice, only : my_order_parameters
use m_Hamiltonian_variables, only : coeff_ham_inter_spec
type(coeff_ham_inter_spec), target, public, protected :: zeeman

private
public :: get_ham_zeeman

contains

subroutine get_ham_zeeman(fname,dim_ham,Ms)
use m_io_files_utils
use m_vector, only : norm
use m_io_utils
use m_convert
use m_constants, only : mu_B,mu_0
implicit none
real(kind=8), intent(in) :: Ms
integer, intent(in) :: dim_ham
character(len=*), intent(in) ::fname
! internal
integer :: io_param
real(kind=8) :: h_ext(3)
! anisotropy
integer :: x_start,x_end
! electric field
integer :: y_start,y_end
integer :: i
character(len=50) :: form

h_ext=0.0d0
zeeman%name='zeeman'
zeeman%order=2

io_param=open_file_read(fname)
call get_parameter(io_param,fname,'c_zeeman',zeeman%c_ham)
! count the magnetic field if present
call get_parameter(io_param,fname,'H_ext',3,h_ext)

call close_file(fname,io_param)

if (norm(H_ext).ge.1.0d-8) then
   zeeman%i_exist=.true.
  else
   return
endif

allocate(zeeman%ham(1))
allocate(zeeman%ham(1)%H(dim_ham,dim_ham))

zeeman%ham(1)%H=0.0d0


call get_borders('magnetic',x_start,x_end,'Bfield',y_start,y_end,my_order_parameters)

! get the zeeman Hamiltonian
call get_symmetric_Op(zeeman%ham(1)%H,mu_0*mu_B*Ms*zeeman%c_ham/2.0d0,x_start,y_start,x_end,y_end)

form=convert('(',dim_ham,'(f12.8,2x))')
write(6,'(a)') ''
write(6,'(a)') 'Zeeman tensor'
do i=1,dim_ham
    write(6,form) zeeman%ham(1)%H(:,i)
enddo
write(6,'(a)') ''

end subroutine get_ham_zeeman


end module m_zeeman
