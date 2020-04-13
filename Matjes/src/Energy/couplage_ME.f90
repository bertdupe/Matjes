module m_couplage_ME
use m_symmetry_operators
use m_lattice, only : my_order_parameters
use m_Hamiltonian_variables, only : coeff_ham_inter_spec
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
real(kind=8), allocatable :: ME_local_sym(:),ME_local_antisym(:),ham_DMI_local(:,:,:)
! magnetization
integer :: x_start,x_end
! electric field
integer :: y_start,y_end
! slope
integer :: i

neighbor_ME_sym=0
neighbor_ME_antisym=0
ME%name='magnetoelectric'

io_param=open_file_read(fname)
call get_parameter(io_param,fname,'c_ME',ME%c_ham)
! count the ME coefficients if present
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
endif

N_ME=max(neighbor_ME_sym,neighbor_ME_antisym)

allocate(ME%ham(N_ME))
do i=1,N_ME
   allocate(ME%ham(i)%H(dim_ham,dim_ham))
   ME%ham(i)%H=0.0d0
enddo

call get_borders('magnetic',x_start,x_end,'Efield',y_start,y_end,my_order_parameters)

! get the symmetric ME effect
do i=1,n_ME
   call get_symmetric_Op(ME%ham(i)%H,ME_local_sym(i)*ME%c_ham/2.0d0,x_start,y_start,x_end,y_end)
enddo

! get the antisymmetric ME effect
allocate(ham_DMI_local(3,3,neighbor_ME_antisym))
ham_DMI_local=0.0
do i=1,neighbor_ME_antisym
    ham_DMI_local(2,1,i)=-ME_local_antisym(i)/2.0d0*ME%c_ham
    ham_DMI_local(3,1,i)=ME_local_antisym(i)/2.0d0*ME%c_ham
    ham_DMI_local(3,2,i)=-ME_local_antisym(i)/2.0d0*ME%c_ham

    ham_DMI_local(1,2,i)=-ham_DMI_local(2,1,i)
    ham_DMI_local(2,3,i)=-ham_DMI_local(3,2,i)
    ham_DMI_local(1,3,i)=-ham_DMI_local(3,1,i)
enddo

do i=1,neighbor_ME_antisym
    call get_Op_in_Op(ME%ham(i)%H,ham_DMI_local(:,:,i),x_start,x_end,y_start,y_end)
    call get_Op_in_Op(ME%ham(i)%H,-ham_DMI_local(:,:,i),y_start,y_end,x_start,x_end)
enddo

end subroutine get_ham_ME

end module m_couplage_ME
