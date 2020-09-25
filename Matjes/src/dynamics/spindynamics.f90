module m_spindynamics
contains
subroutine spindynamics(mag_lattice,mag_motif,io_simu,ext_param,Ham)
use m_basic_types, only : vec_point
use m_derived_types, only : cell,io_parameter,simulation_parameters,point_shell_Operator
use m_type_lattice, only : lattice
use m_modes_variables, only : point_shell_mode
use m_torques, only : get_torques
use m_lattice, only : my_order_parameters
use m_eval_BTeff
use m_measure_temp
use m_topo_commons
use m_update_time
use m_vector, only : cross,norm,norm_cross
use m_randist
use m_constants, only : pi,k_b,hbar
use m_eval_Beff
use m_write_spin
use m_energyfield, only : get_Energy_distrib
use m_createspinfile
!use m_local_energy
use m_dyna_utils
use m_user_info
use m_excitations
use m_operator_pointer_utils
use m_solver_commun
use m_topo_sd
use m_derivative
use m_forces
use m_fftw, only : calculate_fft
use m_plot_FFT
use m_dipolar_field, only : prepare_FFT_dipole,calculate_FFT_modes
use m_solver_order
use m_io_files_utils
use m_tracker
use m_print_Beff
use omp_lib
use m_precision
use m_Htype_gen
implicit none
! input
type(lattice), intent(inout) :: mag_lattice
type(cell), intent(in) :: mag_motif
type(io_parameter), intent(in) :: io_simu
type(simulation_parameters), intent(in) :: ext_param
class(t_H), intent(in) :: Ham
! internal
logical :: gra_log,io_stochafield
integer :: i,j,gra_freq,i_loop,input_excitations
! lattices that are used during the calculations
real(kind=8),allocatable,dimension(:,:,:,:,:,:),target :: spinafter
real(8),pointer     ::  spinafter_1(:)
real(kind=8),allocatable,dimension(:,:) :: D_T,Bini,BT
real(kind=8),allocatable,dimension(:,:,:) :: D_mode
real(kind=8),allocatable,dimension(:) :: D_mode_int
type(lattice)                         :: lat_1,lat_2
! pointers specific for the modes
type(vec_point),allocatable,dimension(:) :: mode_temp,mode_Efield,mode_Hfield,mode_magnetic,mode_disp
type(vec_point),target,allocatable,dimension(:) :: D_T_mag,B_mag,BT_mag
type(vec_point),target,allocatable,dimension(:,:) :: D_mode_mag,mode_excitation_field,lattice_ini_excitation_field
type(vec_point),target,allocatable,dimension(:) :: D_T_disp,B_disp,BT_disp
type(vec_point),target,allocatable,dimension(:,:) :: D_mode_disp
! dummys
real(kind=8) :: qeuler,q_plus,q_moins,vortex(3),Mdy(3),Edy,check1,check2,Eold,check3,Et,dt
real(kind=8) :: Mx,My,Mz,vx,vy,vz,check(2),test_torque,Einitial,ave_torque
real(kind=8) :: dumy(5),security(2)
real(kind=8) :: timestep_int,real_time,h_int(3),damping,E_int(3)
real(kind=8) :: kt,ktini,ktfin
real(kind=8) :: time
integer :: iomp,shape_lattice(4),shape_spin(4),N_cell,N_loop,duration,Efreq
!integer :: io_test
!! switch that controls the presence of magnetism, electric fields and magnetic fields
logical :: i_magnetic,i_temperature,i_mode,i_Efield,i_Hfield,i_excitation,i_displacement
! dumy
logical :: said_it_once,gra_topo

time=0.0d0
input_excitations=0

OPEN(7,FILE='EM.dat',action='write',status='replace',form='formatted')
      Write(7,'(20(a,2x))') '# 1:step','2:real_time','3:E_av','4:M', &
     &  '5:Mx','6:My','7:Mz','8:vorticity','9:vx', &
     &  '10:vy','11:vz','12:qeuler','13:q+','14:q-','15:T=', &
     &  '16:Tfin=','17:Ek=','18:Hx','19:Hy=','20:Hz='

! check the convergence
open(8,FILE='convergence.dat',action='write',form='formatted')

! prepare the matrices for integration

call rw_dyna(timestep_int,damping,Efreq,duration)
shape_lattice=shape(mag_lattice%ordpar%l_modes)
N_cell=product(shape_lattice)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Select the propagators and the integrators
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call select_propagator(ext_param%ktini%value,N_loop)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Create copies of lattice with order-parameter for intermediary states
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Call mag_lattice%copy(lat_1) 
Call mag_lattice%copy(lat_2) 

!shape_spin=shape_lattice
!spinafter=0.0d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! associate pointer for the topological charge, vorticity calculations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call user_info(6,time,'topological operators',.false.)

call get_size_Q_operator(mag_lattice)
call associate_Q_operator(lat_1%ordpar%all_l_modes,mag_lattice%boundary,shape(mag_lattice%ordpar%l_modes))

call user_info(6,time,'done',.true.)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! allocate the element of integrations and associate the pointers to them
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
allocate(D_mode(mag_lattice%dim_mode,N_cell,N_loop),D_T(mag_lattice%dim_mode,N_cell),D_mode_int(mag_lattice%dim_mode))
D_mode=0.0d0
D_T=0.0d0
D_mode_int=0.0d0

allocate(Bini(mag_lattice%dim_mode,N_cell),BT(mag_lattice%dim_mode,N_cell))
Bini=0.0d0
BT=0.0d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! associate pointers only for the magnetization or local modes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
i_magnetic=.false.
i_temperature=.false.
i_mode=.false.
i_Efield=.false.
i_Hfield=.false.

! magnetization
do i=1,size(my_order_parameters)
  if ('magnetic'.eq.trim(my_order_parameters(i)%name)) then
   allocate(mode_magnetic(N_cell),D_mode_mag(N_cell,N_loop),D_T_mag(N_cell),B_mag(N_cell),BT_mag(N_cell))
   call dissociate(mode_magnetic,N_cell)
   call associate_pointer(mode_magnetic,lat_1%ordpar%all_l_modes,'magnetic',i_magnetic)

   do j=1,N_loop
      call dissociate(D_mode_mag(:,j),N_cell)
      call associate_pointer(D_mode_mag(:,j),D_mode(:,:,j),'magnetic',i_magnetic)
   enddo

   call dissociate(D_T_mag,N_cell)
   call associate_pointer(D_T_mag,D_T,'magnetic',i_magnetic)

   call dissociate(B_mag,N_cell)
   call associate_pointer(B_mag,Bini,'magnetic',i_magnetic)

   call dissociate(BT_mag,N_cell)
   call associate_pointer(BT_mag,BT,'magnetic',i_magnetic)

   exit
  endif
enddo

! temperature
do i=1,size(my_order_parameters)
  if ('temperature'.eq.trim(my_order_parameters(i)%name)) then
   allocate(mode_temp(N_cell))
   call dissociate(mode_temp,N_cell)
   call associate_pointer(mode_temp,lat_1%ordpar%all_l_modes,'temperature',i_temperature)

   exit
  endif
enddo

! magnetic field
do i=1,size(my_order_parameters)
  if ('Bfield'.eq.trim(my_order_parameters(i)%name)) then
   allocate(mode_Hfield(N_cell))
   call dissociate(mode_Hfield,N_cell)
   call associate_pointer(mode_Hfield,lat_1%ordpar%all_l_modes,'Bfield',i_Hfield)

   exit
  endif
enddo

! Electric field
do i=1,size(my_order_parameters)
  if ('Efield'.eq.trim(my_order_parameters(i)%name)) then
   allocate(mode_Efield(N_cell))
   call dissociate(mode_Efield,N_cell)
   call associate_pointer(mode_Efield,lat_1%ordpar%all_l_modes,'Efield',i_Efield)

   exit
  endif
enddo

! atomic displacements
do i=1,size(my_order_parameters)
  if ('displacement'.eq.trim(my_order_parameters(i)%name)) then
   allocate(mode_disp(N_cell),D_mode_disp(N_cell,N_loop),D_T_disp(N_cell),B_disp(N_cell),BT_disp(N_cell))
   call associate_pointer(mode_disp,lat_1%ordpar%all_l_modes,'displacement',i_displacement)

   do j=1,N_loop
     call dissociate(D_mode_disp(:,j),N_cell)
     call associate_pointer(D_mode_disp(:,j),D_mode(:,:,j),'displacement',i_displacement)
   enddo

   call dissociate(D_T_disp,N_cell)
   call associate_pointer(D_T_disp,D_T,'displacement',i_displacement)

   call dissociate(B_disp,N_cell)
   call associate_pointer(B_disp,Bini,'displacement',i_displacement)

   call dissociate(BT_disp,N_cell)
   call associate_pointer(BT_disp,BT,'displacement',i_displacement)

   exit
  endif
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Prepare the calculation of the energy and the effective field
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call get_B_matrix(mag_lattice%dim_mode)
!call set_E_matrix(mag_lattice%dim_mode)

#ifdef CPPEIGEN_SPARSE
Call set_large_B(mag_lattice%dim_mode)
spinafter_1(1:dimB)=>spinafter(:,:,:,:,:,1)
#elif defined __sparse_mkl__
Call set_B_sparse(mag_lattice%dim_mode)
spinafter_1(1:dimB)=>spinafter(:,:,:,:,:,1)
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! start the simulation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

gra_log=io_simu%io_Xstruct
io_stochafield=io_simu%io_Tfield
gra_freq=io_simu%io_frequency
gra_topo=io_simu%io_topo
ktini=ext_param%ktini%value
ktfin=ext_param%ktfin%value
kt=ktini
Eold=100.0d0
real_time=0.0d0
Einitial=0.0d0
h_int=ext_param%H_ext%value
E_int=ext_param%E_ext%value
said_it_once=.False.
security=0.0d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! initialize the simulation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Call mag_lattice%copy_val_to(lat_1)

Call Ham%eval_all(Edy,mag_lattice)

write(6,'(a,2x,E20.12E3)') 'Initial Total Energy (eV)',Edy/real(N_cell)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! prepare the derivation of the lattice
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (io_simu%io_Force) call get_derivative(mode_magnetic,mag_lattice)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! prepare the dipole dipole FFT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call prepare_FFT_dipole(N_cell)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! part of the excitations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call get_excitations('input',i_excitation,input_excitations)
! allocate the field on which the excitation occurs
if (i_excitation) then
   allocate(mode_excitation_field(N_cell,input_excitations),lattice_ini_excitation_field(N_cell,input_excitations))
   do iomp=1,input_excitations
     call associate_excitation(iomp,mode_excitation_field(:,iomp),lat_1%ordpar%all_l_modes,my_order_parameters)
     call associate_excitation(iomp,lattice_ini_excitation_field(:,iomp),mag_lattice%ordpar%all_l_modes,my_order_parameters)
   enddo
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! do we use the update timestep
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call init_update_time('input')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! initialize the different torques
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call get_torques('input')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! check if a magnetic texture should be tracked
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (io_simu%io_tracker) call init_tracking(mag_lattice)

!!$OMP parallel
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! beginning of the
do j=1,duration
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   call init_temp_measure(check,check1,check2,check3)

   call truncate(lat_1,mag_lattice%dim_mode)
   qeuler=0.0d0
   q_plus=0.0d0
   q_moins=0.0d0
   vx=0.0d0
   vy=0.0d0
   vz=0.0d0
   Mx=0.0d0
   My=0.0d0
   Mz=0.0d0
   Edy=0.0d0
   Mdy=0.0d0
   vortex=0.0d0
   test_torque=0.0d0
   ave_torque=0.0d0
   Mdy=0.0d0
!!$OMP do private(iomp) schedule(auto)
   do iomp=1,N_cell
     Bini(:,iomp)=0.0d0
     D_mode(:,iomp,:)=0.0d0
     BT(:,iomp)=0.0d0
   enddo
!!$OMP end do
   dt=timestep_int

   call update_ext_EM_fields(real_time,check)

   call calculate_FFT_modes(j)

   test_torque=0.0d0



!
! loop over the integration order
!
!!$OMP do ordered private(i_loop,iomp) schedule(auto)
  do i_loop=1,N_loop

!
! loop that get all the fields
!
    do iomp=1,N_cell
      if (i_excitation) then
        call update_EMT_of_r(iomp,mode_excitation_field)
        call update_EMT_of_r(iomp,lattice_ini_excitation_field)
      endif
    enddo 

#ifdef CPPEIGEN_SPARSE
    Call energy_B(Bini,size(Bini,1),size(Bini,2),spinafter_1)
#elif defined __sparse_mkl__
    Call B_sparse(Bini,size(Bini,1),size(Bini,2),spinafter_1) 
#else
    Bini=0.0d0
    do iomp=1,N_cell
      call calculate_Beff(Bini(:,iomp),iomp,lat_1%ordpar%all_l_modes)
    enddo 
#endif

!
! Be carefull the sqrt(dt) is not included in BT_mag(iomp),D_T_mag(iomp) at this point. It is included only during the integration
!
    do iomp=1,N_cell
      if (i_temperature) call get_temperature_field(mode_temp(iomp)%w(1),damping,mode_magnetic(iomp)%w,BT_mag(iomp)%w,D_T_mag(iomp)%w,size(mode_magnetic(iomp)%w))
      if (i_magnetic) D_mode_mag(iomp,i_loop)%w=get_propagator_field(B_mag(iomp)%w,damping,mode_magnetic(iomp)%w,size(mode_magnetic(iomp)%w))
    enddo

!
! check the Butcher's table for the new dt
! if dt is 0 in the Butcher's table do not carry out integration
    dt=get_dt_mode(timestep_int,i_loop)
!
! loop that carry out the integration
!

    do iomp=1,N_cell

       call get_D_mode(D_mode(:,iomp,:),i_loop,N_loop,D_mode_int)

       lat_2%ordpar%all_l_modes(iomp)%w=get_integrator_field(mag_lattice%ordpar%all_l_modes(iomp)%w,D_mode_int,D_T(:,iomp),dt,mag_lattice%dim_mode)

    enddo

!!!!!!!!!!!!!!! copy the final configuration from spinafter(:,2) in spinafter(:,1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


   Call lat_2%copy_val_to(lat_1)

  enddo
!!$OMP end do
!!$OMP single
!!!!!!!!!!!!!!! copy the final configuration in my_lattice
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Call lat_2%copy_val_to(mag_lattice)

call truncate(mag_lattice,mag_lattice%dim_mode)

!
!!!!!! Measure the temperature if the users wish
!
!!!!$omp do private(iomp) reduction(+:check1,check2) reduction(max:test_torque) schedule(auto)
if (i_temperature) then
  do iomp=1,N_cell

    call update_temp_measure(check1,check2,mode_magnetic(iomp)%w,B_mag(iomp)%w)
    if (norm_cross(mode_magnetic(iomp)%w,B_mag(iomp)%w,1,3).gt.test_torque) test_torque=norm_cross(mode_magnetic(iomp)%w,B_mag(iomp)%w,1,3)
  enddo
endif
!!!!$omp end do
check(1)=check(1)+check1
check(2)=check(2)+check2

#ifdef CPP_OPENMP
!!$OMP end parallel
#endif


if (j.eq.1) check3=test_torque

! calculate energy

!#ifdef CPPEIGEN_SPARSE
!Call energy_H(Edy,mag_lattice%dim_mode)
!#else
Call Ham%eval_all(Edy,mag_lattice)
!#endif
do iomp=1,N_cell

    Mdy=Mdy+mode_magnetic(iomp)%w

    dumy=get_charge(iomp)

    q_plus=q_plus+dumy(1)/pi(4.0d0)
    q_moins=q_moins+dumy(2)/pi(4.0d0)

    vx=vx+dumy(3)
    vy=vy+dumy(4)
    vz=vz+dumy(5)

enddo

!!!$omp end do
vortex=(/vx,vy,vz/)/3.0d0/sqrt(3.0d0)
Edy=Edy/real(N_cell)
Mdy=Mdy/real(N_cell)

if (dabs(check(2)).gt.1.0d-8) call get_temp(security,check,kt)

!
! update pattern recognition
!

if (io_simu%io_tracker) then
!  call update_tracking(j)
  if (mod(j-1,gra_freq).eq.0) call plot_tracking(j/gra_freq,lat_1)
endif


!!!!!!!!!!!!!!! write local effective field to file here
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if ((io_simu%io_Beff).and.(mod(j-1,gra_freq).eq.0)) call print_Beff(j/gra_freq,Bini)

if (mod(j-1,Efreq).eq.0) Write(7,'(I6,18(E20.12E3,2x),E20.12E3)') j,real_time,Edy, &
     &   norm(Mdy),Mdy,norm(vortex),vortex,q_plus+q_moins,q_plus,q_moins, &
     &   kT/k_B,(security(i),i=1,2),H_int

if ((io_simu%io_Energy_Distrib).and.((mod(j-1,gra_freq).eq.0))) then
    call get_Energy_distrib(j/gra_freq,mag_lattice%ordpar%all_l_modes)
endif

if ((gra_log).and.(mod(j-1,gra_freq).eq.0)) then
    call CreateSpinFile(j/gra_freq,mag_lattice%ordpar%all_l_modes)
    call WriteSpinAndCorrFile(j/gra_freq,mag_lattice%ordpar%all_l_modes,'SpinSTM_')
    write(6,'(a,3x,I10)') 'wrote Spin configuration and povray file number',j/gra_freq
    write(6,'(a,3x,f14.6,3x,a,3x,I10)') 'real time in ps',real_time/1000.0d0,'iteration',j
endif

if ((io_stochafield).and.(mod(j-1,gra_freq).eq.0)) then
    call WriteSpinAndCorrFile(j/gra_freq,BT,'Stocha-field_')
    write(6,'(a,I10)')'wrote Spin configuration and povray file number',j/gra_freq
endif

if ((gra_topo).and.(mod(j-1,gra_freq).eq.0)) Call get_charge_map(j/gra_freq)

if ((io_simu%io_Force).and.(mod(j-1,gra_freq).eq.0)) call forces(j/gra_freq,lat_1%ordpar%all_l_modes,mag_lattice%dim_mode,mag_lattice%areal)

if ((io_simu%io_fft_Xstruct).and.(mod(j-1,gra_freq).eq.0)) call plot_fft(mag_lattice%ordpar%all_l_modes,-1.0d0,mag_lattice%areal,mag_lattice%dim_lat,mag_lattice%boundary,mag_lattice%dim_mode,j/gra_freq)

! security in case of energy increase in SD and check for convergence
if (((damping*(Edy-Eold).gt.1.0d-10).or.(damping*(Edy-Einitial).gt.1.0d-10)).and.(kt.lt.1.0d-10).and.(.not.said_it_once)) then
#ifdef CPP_MPI
    write(6,'(a)') 'WARNING: the total energy or torque is increasing for non zero damping'
    write(6,'(a)') 'this is not allowed by theory'
    write(6,'(a)') 'please reduce the time step'
#else
    write(6,'(a)') 'WARNING: the total energy or torque is increasing for non zero damping'
    write(6,'(a)') 'this is not allowed by theory'
    write(6,'(a)') 'please reduce the time step'
#endif
    said_it_once=.True.
endif

if (mod(j-1,Efreq).eq.0) write(8,'(I10,3x,3(E20.12E3,3x))') j,Edy,test_torque,ave_torque

check3=test_torque
Eold=Edy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! update timestep
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call update_time(timestep_int,Bini,BT,damping)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! reinitialize T variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
check(1)=0.0d0
check(2)=0.0d0


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!! end of a timestep
real_time=real_time+timestep_int !increment time counter
enddo 

#ifdef CPP_OPENMP
!!$OMP end parallel
#endif
!!!!!!!!!!!!!!! end of a timestep
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef CPP_MPI
      if (irank.eq.0) then
#endif
      close(7)
      close(8)
#ifdef CPP_MPI
      endif
#endif

if ((dabs(check(2)).gt.1.0d-8).and.(kt/k_B.gt.1.0d-5)) then
    write(6,'(a,2x,f16.6)') 'Final Temp (K)', check(1)/check(2)/2.0d0/k_B
    write(6,'(a,2x,f14.7)') 'Kinetic energy (meV)', (check(1)/check(2)/2.0d0-kT)/k_B*1000.0d0
else
    write(6,'(a)') 'the temperature measurement is not possible'
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  copy the last spin lattice into the mag_lattice
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



call kill_B_matrix()
!call kill_E_matrix()


end subroutine spindynamics
end module
