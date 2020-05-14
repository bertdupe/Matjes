module m_minimize

real(kind=8) :: masse=1.0d0
real(kind=8) :: dt=0.1d0
integer :: N_minimization=1000
real(kind=8) :: conv_torque=1.0d-6

interface minimize
  module procedure minimize_2Dlattice,minimize_lattice
end interface

private
public :: minimize

contains

!
!
! Small subroutine that reads the actual lattice
! and send it into minimize_2Dlattice
!
!
subroutine minimize_lattice(my_lattice,io_simu)
use m_derived_types, only : io_parameter,lattice
implicit none
type(io_parameter), intent(in) :: io_simu
type(lattice), intent(inout) :: my_lattice
! internal
real(kind=8), allocatable, dimension(:,:) :: internal_matrix
integer :: shape_lattice(4),dim_mode,N_cell,i,j,k,l,n

shape_lattice=shape(my_lattice%l_modes)
dim_mode=my_lattice%dim_mode
N_cell=product(shape_lattice)

allocate(internal_matrix(dim_mode,N_cell))
internal_matrix=0.0d0

n=0
do l=1,shape_lattice(4)
 do k=1,shape_lattice(3)
   do j=1,shape_lattice(2)
     do i=1,shape_lattice(1)

     n=n+1
     internal_matrix(:,n)=my_lattice%l_modes(i,j,k,l)%w

     enddo
   enddo
 enddo
enddo

call minimize_2Dlattice(internal_matrix,io_simu)

n=0
do l=1,shape_lattice(4)
 do k=1,shape_lattice(3)
   do j=1,shape_lattice(2)
     do i=1,shape_lattice(1)

     n=n+1
     my_lattice%l_modes(i,j,k,l)%w=internal_matrix(:,n)

     enddo
   enddo
 enddo
enddo

end subroutine


!
!
! Minimization routine that it actually used.
! The interface is used to put the data into the good format
!
!
subroutine minimize_2Dlattice(my_lattice,io_simu)
use m_derived_types, only : io_parameter
use m_basic_types, only : vec_point
use m_write_spin
use m_createspinfile
use m_solver, only : minimization
use m_vector, only : norm_cross,norm, calculate_damping
use m_dyna_utils, only : copy_lattice
use m_eval_Beff
use m_local_energy
use m_lattice, only : my_order_parameters
use m_operator_pointer_utils
use omp_lib
implicit none
type(io_parameter), intent(in) :: io_simu
real(kind=8), intent(inout) :: my_lattice(:,:)
! dummy variable
real(kind=8),allocatable, dimension(:,:) :: velocity,predicator,force
real(kind=8),allocatable, dimension(:) :: F_eff,V_eff,F_temp
type(vec_point),allocatable,dimension(:) :: all_mode,mode_magnetic
! internal
real(kind=8) :: dumy,force_norm,Energy,vmax,vtest,Eint,test_torque,max_torque
! the computation time
integer :: i_min,gra_freq,i,shape_lattice(2)
logical :: gra_log,i_magnetic
integer :: iomp,dim_mode,N_cell
#ifdef CPP_OPENMP
integer :: nthreads,ithread
#endif

call init_variables()

gra_freq=io_simu%io_frequency
gra_log=io_simu%io_Xstruct
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

force_norm=0.0d0
test_torque=0.0d0
max_torque=10.0d0
vmax=0.0d0
vtest=0.0d0
shape_lattice=shape(my_lattice)
N_cell=shape_lattice(2)
dim_mode=shape_lattice(1)

allocate(all_mode(N_cell))
call associate_pointer(all_mode,my_lattice)

allocate(velocity(dim_mode,N_cell),predicator(dim_mode,N_cell),force(dim_mode,N_cell))
allocate(F_eff(dim_mode),V_eff(dim_mode),F_temp(dim_mode))
velocity=0.0d0
force=0.0d0
predicator=0.0d0
F_eff=0.0d0
V_eff=0.0d0
F_temp=0.0d0

! magnetization
do i=1,size(my_order_parameters)
  if ('magnetic'.eq.trim(my_order_parameters(i)%name)) then
   allocate(mode_magnetic(N_cell))
   call dissociate(mode_magnetic,N_cell)
   call associate_pointer(mode_magnetic,all_mode,'magnetic',i_magnetic)

   exit
  endif
enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Prepare the calculation of the energy and the effective field
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call get_B_matrix(dim_mode)
call get_E_matrix(dim_mode)


do iomp=1,N_cell

   call calculate_Beff(F_eff,iomp,all_mode)

   force(:,iomp)=calculate_damping(all_mode(iomp)%w,F_eff)
   call minimization(all_mode(iomp)%w,force(:,iomp),predicator(:,iomp),dt**2,masse*2.0d0)

   test_torque=norm_cross(predicator(:,iomp),force(:,iomp),1,3)

enddo
call copy_lattice(predicator,all_mode)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! end of initialization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do i_min=1,N_minimization

  max_torque=0.0d0
  dumy=0.0d0
  force_norm=0.0d0
  vmax=0.0d0

  do iomp=1,N_cell

    call calculate_Beff(F_eff,iomp,all_mode)
    F_temp=calculate_damping(all_mode(iomp)%w,F_eff)

    call minimization(velocity(:,iomp),(force(:,iomp)+F_temp)/2.0d0,V_eff,dt,masse)

    F_eff=F_temp
    force(:,iomp)=F_eff
    velocity(:,iomp)=V_eff

    dumy=dumy+dot_product(V_eff,F_eff)
    force_norm=force_norm+norm(F_eff)**2

  enddo

  if (abs(dumy).gt.1.0d-8) then
    do iomp=1,N_cell

      velocity(:,iomp)=dumy*force(:,iomp)/force_norm

   enddo
  else
     velocity=0.0d0
  endif

  do iomp=1,N_cell

    call minimization(all_mode(iomp)%w,velocity(:,iomp),force(:,iomp),predicator(:,iomp),dt,masse)
    test_torque=norm(force(:,iomp))
    vtest=norm(velocity(:,iomp))**2

    if (vtest.gt.vmax) vmax=vtest

    if (test_torque.gt.max_torque) max_torque=test_torque

  enddo

  call copy_lattice(predicator,all_mode)

  Energy=0.0d0
  do iomp=1,N_cell

    call local_energy(Eint,iomp,all_mode)
    Energy=Energy+Eint

  enddo

  write(6,'(/,a,2x,I10)') 'iteration',i_min
  write(6,'(a,2x,f14.11)') 'Energy of the system (eV/unit cell)',Energy/dble(N_cell)
  write(6,'(2(a,2x,f14.11,2x))') 'convergence criteria:',conv_torque,',Measured Torque:',max_torque
  write(6,'(a,2x,f14.11,/)') 'speed of displacements:',vmax

  if ((gra_log).and.(mod(i_min-1,gra_freq).eq.0)) then
    call WriteSpinAndCorrFile((i_min-1)/gra_freq,all_mode,'spin_minimization')
    call CreateSpinFile((i_min-1)/gra_freq,all_mode)
  endif

  if (conv_torque.gt.max_torque) then
    write(6,'(a)') 'minimization converged'
    exit
  endif

enddo ! number of minimization steps


call kill_B_matrix()
call kill_E_matrix()

end subroutine


!!
!
! Initialize the energy minimization
!
!!
subroutine init_variables()
use m_io_files_utils
use m_io_utils
! internal
integer :: io


io=open_file_read('input')

call get_parameter(io,'input','steps',N_minimization)
call get_parameter(io,'input','masse',masse)
call get_parameter(io,'input','timestep',dt)
call get_parameter(io,'input','convergence_criteria',conv_torque)

call close_file('input',io)

!!! check the input
if (masse.eq.0.0d0) then
   write(6,'(a)') 'The mass should be different from 0'
   stop
endif

end subroutine

end module m_minimize
