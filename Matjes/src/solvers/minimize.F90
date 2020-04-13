module m_minimize

real(kind=8) :: masse=1.0d0
real(kind=8) :: dt=0.1d0
integer :: N_minimization=1000
real(kind=8) :: conv_torque=1.0d-6

private

public :: minimize

contains



subroutine minimize(spin,mode_B_column,B_line,mode_E_column,E_line,io_simu)
use m_derived_types, only : point_shell_Operator, io_parameter
use m_modes_variables, only : point_shell_mode
use m_basic_types, only : vec_point
use m_write_spin
use m_createspinfile
use m_solver, only : minimization
use m_vector, only : norm_cross,norm, calculate_damping
use m_dyna_utils, only : copy_lattice
use m_eval_Beff
use m_local_energy
#ifdef CPP_OPENMP
      use omp_lib
#endif
implicit none
type(io_parameter), intent(in) :: io_simu
type(vec_point), intent(inout) :: spin(:)
type(point_shell_Operator), intent(in), dimension(:) :: B_line,E_line
type(point_shell_mode), intent(in), dimension(:) :: mode_B_column,mode_E_column
! dummy variable
real(kind=8),allocatable, dimension(:,:) :: velocity,predicator,force
!real(kind=8) :: velocity(3,shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5))
!real(kind=8) :: predicator(3,shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5))
!real(kind=8) :: force(3,shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5))

! internal
real(kind=8) :: F_eff(3),V_eff(3),dumy,force_norm,Energy,vmax,vtest,F_temp(3),Eint,test_torque,max_torque
! the computation time
integer :: i_min,N_cell,gra_freq
logical :: gra_log
integer :: iomp
#ifdef CPP_OPENMP
integer :: nthreads,ithread
#endif

call init_variables()

gra_freq=io_simu%io_frequency
gra_log=io_simu%io_Xstruct
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

N_cell=size(spin)
allocate(velocity(3,N_cell),predicator(3,N_cell),force(3,N_cell))
velocity=0.0d0
force=0.0d0
predicator=0.0d0
F_eff=0.0d0
V_eff=0.0d0
force_norm=0.0d0
test_torque=0.0d0
max_torque=10.0d0
vmax=0.0d0
vtest=0.0d0

#ifdef CPP_OPENMP
nthreads=omp_get_num_procs()
call omp_set_num_threads(nthreads)

write(6,'(a)') 'OMP parallelization selected'
write(6,'(a,2x,I3,2x,a)') 'I will calculate on',nthreads,'threads'

#endif

do iomp=1,N_cell

   call calculate_Beff(F_eff,mode_B_column(iomp),B_line(iomp),iomp)
   force(:,iomp)=calculate_damping(spin(iomp)%w,F_eff)
   call minimization(spin(iomp)%w,force(:,iomp),predicator(:,iomp),dt**2,masse*2.0d0,3)

   test_torque=norm_cross(predicator(:,iomp),force(:,iomp),1,3)
   if (test_torque.gt.max_torque) max_torque=test_torque

enddo

call copy_lattice(predicator,spin)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! end of initialization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(*,*) N_minimization
do i_min=1,N_minimization

max_torque=0.0d0
dumy=0.0d0
force_norm=0.0d0
vmax=0.0d0

#ifdef CPP_OPENMP
!$OMP parallel default(shared) private(ithread)
ithread=omp_get_thread_num()

!$OMP do private(F_eff,V_eff,iomp) reduction(+:dumy,force_norm) reduction(max:vmax) collapse(4)

#endif

do iomp=1,N_cell

   call calculate_Beff(F_eff,mode_B_column(iomp),B_line(iomp),iomp)
   F_temp=calculate_damping(spin(iomp)%w,F_eff)

   call minimization(velocity(:,iomp),(force(:,iomp)+F_temp)/2.0d0,V_eff,dt,masse,3)

   F_eff=F_temp
   force(:,iomp)=F_eff
   velocity(:,iomp)=V_eff

   dumy=dumy+dot_product(V_eff,F_eff)
   force_norm=force_norm+norm(F_eff)

enddo

#ifdef CPP_OPENMP
!$OMP end do
#endif

if (dumy.gt.0.0d0) then
#ifdef CPP_OPENMP
!$OMP do collapse(1)
#endif
   do iomp=1,N_cell

      velocity(:,iomp)=dumy*force(:,iomp)/force_norm

   enddo
#ifdef CPP_OPENMP
!$OMP end do
#endif
else
   velocity=0.0d0
endif

#ifdef CPP_OPENMP
!$OMP do private(iomp,test_torque) reduction(+:max_torque) collapse(1)
#endif
do iomp=1,N_cell

   call minimization(spin(iomp)%w,velocity(:,iomp),force(:,iomp),predicator(:,iomp),dt,masse,3)
   test_torque=norm(force(:,iomp))
   vtest=norm(velocity(:,iomp))**2

   if (vtest.gt.vmax) vmax=vtest

   if (test_torque.gt.max_torque) max_torque=test_torque

enddo

#ifdef CPP_OPENMP
!$OMP end do
#endif

call copy_lattice(predicator,spin)
Energy=0.0d0

#ifdef CPP_OPENMP
!$OMP do private(iomp) reduction(+:Energy) collapse(1)
#endif

do iomp=1,N_cell

   call local_energy(Eint,iomp,mode_E_column(iomp),E_line(iomp))
   Energy=Energy+Eint

enddo

#ifdef CPP_OPENMP
!$OMP end do
!$OMP end parallel
#endif

       

if ((gra_log).and.(mod(i_min-1,gra_freq).eq.0)) then
   write(6,'(/,a,2x,I10)') 'iteration',i_min
   write(6,'(a,2x,f14.11)') 'Energy of the system (eV/unit cell)',Energy/dble(N_cell)
   write(6,'(2(a,2x,f14.11,2x))') 'convergence criteria:',conv_torque,',Measured Torque:',max_torque
   write(6,'(a,2x,f14.11,/)') 'speed of displacements:',vmax
   call WriteSpinAndCorrFile(i_min/gra_freq,spin,'spin_minimization')
   call CreateSpinFile(i_min/gra_freq,spin)
endif

if (conv_torque.gt.max_torque) then
   write(6,'(a)') 'minimization converged'
   exit
endif

enddo ! number of minimization steps

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
