module m_exchange_heisenberg
use m_symmetry_operators
use m_lattice, only : my_order_parameters
use m_Hamiltonian_variables, only : coeff_ham_inter_spec
type(coeff_ham_inter_spec), target, public, protected :: exchange

private
public :: get_ham_exchange

contains

subroutine get_ham_exchange(fname,dim_ham)
use m_io_files_utils
use m_io_utils
use m_convert
implicit none
integer, intent(in) :: dim_ham
character(len=*), intent(in) ::fname
! internal
integer :: neighbor_exch_sym,io_param,neighbor_exch_antisym,N_exch
real(kind=8), allocatable :: exch_local_sym(:),exch_local_antisym(:),ham_DMI_local(:,:,:),c_DMI
! magnetization
integer :: x_start,x_end
! electric field
integer :: y_start,y_end
! slope
integer :: i,j
character(len=50) :: form

c_DMI=-1.0d0
neighbor_exch_sym=0
neighbor_exch_antisym=0
exchange%name='exchange'
exchange%N_shell=-1
exchange%order=2
exch_local_sym=0.0d0

io_param=open_file_read(fname)
call get_parameter(io_param,fname,'c_Jij',exchange%c_ham)
! count the exchange coefficients if present
! count the number of anisotropy coefficients if present
call get_parameter(io_param,fname,'N_exchange',exchange%N_shell)
if (exchange%N_shell.eq.-1) then
   neighbor_exch_sym=count_variables(io_param,'J_',fname)
else
   neighbor_exch_sym=exchange%N_shell
endif


if (neighbor_exch_sym.ne.0) then
   allocate(exch_local_sym(neighbor_exch_sym))
   exch_local_sym=0.0d0

   call get_coeff(io_param,fname,'J_',exch_local_sym)
   neighbor_exch_sym=number_nonzero_coeff(exch_local_sym,'exchange')
endif
if (neighbor_exch_sym.ne.0) exchange%i_exist=.true.

neighbor_exch_antisym=count_variables(io_param,'DMI_',fname)
if (neighbor_exch_antisym.ne.0) then
   allocate(exch_local_antisym(neighbor_exch_antisym))
   exch_local_antisym=0.0d0

   call get_coeff(io_param,fname,'DMI_',exch_local_antisym)
   neighbor_exch_antisym=number_nonzero_coeff(exch_local_antisym,'DMI')
endif

if ((neighbor_exch_antisym.eq.0).and.(neighbor_exch_sym.eq.0)) then
   return
  else
   exchange%i_exist=.true.
endif

N_exch=max(neighbor_exch_sym,neighbor_exch_antisym)

!
! if no exchange is present just return
!
if (N_exch.eq.0) return


allocate(exchange%ham(N_exch))
do i=1,N_exch
   allocate(exchange%ham(i)%H(dim_ham,dim_ham))
   exchange%ham(i)%H=0.0d0
enddo

call get_borders('magnetic',x_start,x_end,'magnetic',y_start,y_end,my_order_parameters)

! get the symmetric exchange
if (neighbor_exch_sym.ne.0) then
  call get_parameter(io_param,fname,'c_DM',c_DMI)
  do i=1,n_exch
     call get_diagonal_Op(exchange%ham(i)%H,exch_local_sym(i),exchange%c_ham,x_start,x_end)
  enddo
endif

if (neighbor_exch_antisym.ne.0) then

  allocate(ham_DMI_local(3,3,neighbor_exch_antisym))
  ham_DMI_local=0.0d0
  do i=1,neighbor_exch_antisym
    ham_DMI_local(2,1,i)=exch_local_antisym(i)*c_DMI
    ham_DMI_local(3,1,i)=-exch_local_antisym(i)*c_DMI
    ham_DMI_local(3,2,i)=exch_local_antisym(i)*c_DMI

    ham_DMI_local(1,2,i)=-ham_DMI_local(2,1,i)
    ham_DMI_local(2,3,i)=-ham_DMI_local(3,2,i)
    ham_DMI_local(1,3,i)=-ham_DMI_local(3,1,i)
  enddo

  do i=1,neighbor_exch_antisym
    call get_Op_in_Op(exchange%ham(i)%H,ham_DMI_local(:,:,i),x_start,x_end,y_start,y_end)
  enddo

endif

call close_file(fname,io_param)

form=convert('(',dim_ham,'(f12.8,2x))')
write(6,'(a)') ''
write(6,'(a)') 'Exchange tensor of order 2'
do i=1,N_exch
  write(6,'(a,I3)') 'Shell', i
  do j=1,dim_ham
    write(6,form) exchange%ham(i)%H(:,j)
  enddo
enddo
write(6,'(a)') ''

end subroutine get_ham_exchange

end module m_exchange_heisenberg
