module m_anisotropy_heisenberg
use m_symmetry_operators
use m_lattice, only : my_order_parameters
use m_Hamiltonian_variables, only : coeff_ham_inter_spec
type(coeff_ham_inter_spec), target, public, protected :: anisotropy

private
public :: get_ham_anisotropy

contains

subroutine get_ham_anisotropy(fname,dim_ham)
use m_io_files_utils
use m_io_utils
use m_convert
implicit none
integer, intent(in) :: dim_ham
character(len=*), intent(in) ::fname
! internal
integer :: neighbor_ani,io_param
real(kind=8), allocatable :: ani_local_sym(:,:),ham_local(:,:,:)
! anisotropy
integer :: x_start,x_end
! electric field
integer :: y_start,y_end
! slope
integer :: i
character(len=50) :: form

neighbor_ani=0
anisotropy%name='anisotropy'
anisotropy%c_ham=1.0d0
anisotropy%N_shell=-1
anisotropy%order=2

io_param=open_file_read(fname)
call get_parameter(io_param,fname,'c_ani',anisotropy%c_ham)
! count the number of anisotropy coefficients if present
call get_parameter(io_param,fname,'N_ani',anisotropy%N_shell)
if (anisotropy%N_shell.eq.-1) then
   neighbor_ani=count_variables(io_param,'ani_',fname)
else
   neighbor_ani=anisotropy%N_shell
endif


if (neighbor_ani.ne.0) then
   allocate(ani_local_sym(3,neighbor_ani))
   ani_local_sym=0.0d0

   do i=1,neighbor_ani
      call get_parameter(io_param,fname,'ani_',3,ani_local_sym(:,1))
   enddo
   neighbor_ani=number_nonzero_coeff(ani_local_sym,'anisotropy')
endif

call close_file(fname,io_param)

if (neighbor_ani.ne.0) then
   anisotropy%i_exist=.true.
  else
   return
endif

allocate(anisotropy%ham(neighbor_ani))
do i=1,neighbor_ani
   allocate(anisotropy%ham(i)%H(dim_ham,dim_ham))
   anisotropy%ham(i)%H=0.0d0
enddo

call get_borders('magnetic',x_start,x_end,'magnetic',y_start,y_end,my_order_parameters)

! get the symmetric exchange
allocate(ham_local(dim_ham,dim_ham,neighbor_ani))
ham_local=0.0d0
call get_diagonal_Op(ham_local,ani_local_sym,anisotropy%c_ham,x_start,x_end)

do i=1,neighbor_ani
   anisotropy%ham(i)%h=ham_local(:,:,i)
enddo

form=convert('(',dim_ham,'(f12.8,2x))')
write(6,'(a)') ''
write(6,'(a)') 'Anisotropy tensor of order 2'
do i=1,dim_ham
    write(6,form) anisotropy%ham(1)%H(:,i)
enddo
write(6,'(a)') ''

end subroutine get_ham_anisotropy

end module m_anisotropy_heisenberg
