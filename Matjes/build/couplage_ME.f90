module m_couplage_ME
use m_symmetry_operators
use m_lattice, only : my_order_parameters
use m_Hamiltonian_variables, only : coeff_ham_inter_spec
use m_convert
type(coeff_ham_inter_spec), target, public, protected :: ME

private
public :: get_ham_ME,get_number_EM_DMI

contains


integer function get_number_EM_DMI(fname)
use m_io_files_utils
use m_io_utils
implicit none
character(len=*) :: fname
! internal variables
integer :: io

io=open_file_read(fname)
get_number_EM_DMI=count_variables(io,'ME_antisym_',fname)
call close_file(fname,io)

end function get_number_EM_DMI

subroutine get_ham_ME(fname,dim_ham)
use m_io_files_utils
use m_io_utils
implicit none
integer, intent(in) :: dim_ham
character(len=*), intent(in) ::fname
! internal
integer :: neighbor_ME_sym,io_param,neighbor_ME_antisym,N_ME
real(kind=8), allocatable :: ME_local_sym(:),ME_local_antisym(:)
! magnetization
integer :: x_start,x_end
! electric field
integer :: y_start,y_end
! slope
integer :: i,j
character(len=50) :: form

neighbor_ME_sym=0
neighbor_ME_antisym=0
ME%name='magnetoelectric'
ME%order=3

io_param=open_file_read(fname)
call get_parameter(io_param,fname,'c_ME',ME%c_ham)
!
! count the ME coefficients if present
!
neighbor_ME_sym=count_variables(io_param,'ME_sym_',fname)
if (neighbor_ME_sym.ne.0) then
   allocate(ME_local_sym(neighbor_ME_sym))
   ME_local_sym=0.0d0

   call get_coeff(io_param,fname,'ME_sym_',ME_local_sym)
   neighbor_ME_sym=number_nonzero_coeff(ME_local_sym,'ME_sym')
endif

neighbor_ME_antisym=count_variables(io_param,'ME_antisym_',fname)
if (neighbor_ME_antisym.ne.0) then
   allocate(ME_local_antisym(neighbor_ME_antisym))
   ME_local_antisym=0.0d0

   call get_coeff(io_param,fname,'ME_antisym_',ME_local_antisym)
   neighbor_ME_antisym=number_nonzero_coeff(ME_local_antisym,'ME_antisym')
endif

call close_file(fname,io_param)

if ((neighbor_ME_antisym.eq.0).and.(neighbor_ME_sym.eq.0)) then
   return
  else
   ME%i_exist=.true.
   write(6,'(a)') 'WARNING!!!! ME coupling found. You need E_ext different from 0'
endif

N_ME=max(neighbor_ME_sym,neighbor_ME_antisym)

allocate(ME%ham(N_ME))
do i=1,N_ME
   allocate(ME%ham(i)%H(dim_ham,dim_ham**2))
   ME%ham(i)%H=0.0d0
enddo

call get_borders('magnetic',x_start,x_end,'Efield',y_start,y_end,my_order_parameters)

! get the symmetric ME effect
do i=1,n_ME
  ! get diagonal terms
  ME%ham(i)%H(y_start,x_start)=ME_local_sym(i)*ME%c_ham
  ME%ham(i)%H(y_start+1,x_start+dim_ham+1)=ME_local_sym(i)*ME%c_ham
  ME%ham(i)%H(y_start+2,x_start+2*(dim_ham+1))=ME_local_sym(i)*ME%c_ham

  ! get the of diagonal terms
  ME%ham(i)%H(y_start,x_start+1)=-ME_local_antisym(i)*ME%c_ham
  ME%ham(i)%H(y_start+1,x_start+1)=-ME_local_antisym(i)*ME%c_ham

  ME%ham(i)%H(y_start,x_start+2)=ME_local_antisym(i)*ME%c_ham
  ME%ham(i)%H(y_start+2,x_start+2)=ME_local_antisym(i)*ME%c_ham

  ME%ham(i)%H(y_start,x_start+dim_ham)=ME_local_antisym(i)*ME%c_ham
  ME%ham(i)%H(y_start+1,x_start+dim_ham)=ME_local_antisym(i)*ME%c_ham

  ME%ham(i)%H(y_start+1,x_start+dim_ham+2)=-ME_local_antisym(i)*ME%c_ham
  ME%ham(i)%H(y_start+2,x_start+dim_ham+2)=-ME_local_antisym(i)*ME%c_ham

  ME%ham(i)%H(y_start,x_start+2*dim_ham)=-ME_local_antisym(i)*ME%c_ham
  ME%ham(i)%H(y_start+2,x_start+2*dim_ham)=-ME_local_antisym(i)*ME%c_ham

  ME%ham(i)%H(y_start+1,x_start+2*dim_ham+1)=ME_local_antisym(i)*ME%c_ham
  ME%ham(i)%H(y_start+2,x_start+2*dim_ham+1)=ME_local_antisym(i)*ME%c_ham
enddo

form=convert('(',dim_ham,'(f12.8,2x))')

write(6,'(a)') ''
write(6,'(a)') 'Magnetoelectric Hamiltonian of order 3'
do i=1,N_ME
  write(6,'(a,I3)') 'Shell  ',i
  do j=1,dim_ham**2
    write(6,form) ME%ham(i)%H(:,j)
  enddo
  write(6,'(a)') ''
enddo

end subroutine get_ham_ME

end module m_couplage_ME
