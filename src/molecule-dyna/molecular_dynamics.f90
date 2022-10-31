!
!subroutine that carries out molecular dynamics
module m_molecular_dynamics
use m_derived_types, only : lattice,number_different_order_parameters
use m_derived_types, only : io_parameter,simulation_parameters
use m_hamiltonian_collection, only: hamiltonian
use mpi_basic, only: mpi_type
use m_H_public
use m_random_public
use,intrinsic :: iso_fortran_env, only : error_unit, output_unit

implicit none
contains

subroutine molecular_dynamics(my_lattice,io_simu,ext_param,H,comm)
   type(lattice), intent(inout)                 :: my_lattice
   type(io_parameter), intent(inout)            :: io_simu
   type(simulation_parameters), intent(inout)   :: ext_param
   type(hamiltonian),intent(inout)              :: H
   type(mpi_type),intent(in)                    :: comm

   !molecular dynamics so far only implemented without mpi_parallelization
   if(comm%ismas)then
       write(output_unit,'(a)') 'entering into the molecular dynamics routines'
       if(comm%Np>1)then
          write(error_unit,'(//A)') "WARNING, running Molecular dynamics with MPI parallelization which is not implemented."
          write(error_unit,'(A)')   "         Only one thread will to anything."
       endif
       Call molecular_dynamics_run(my_lattice,io_simu,ext_param,H)
   endif
end subroutine


subroutine molecular_dynamics_run(my_lattice,io_simu,ext_param,H)
   use m_derived_types, only : lattice,io_parameter,number_different_order_parameters,simulation_parameters
   use m_energy_output_contribution, only:Eout_contrib_init, Eout_contrib_write
   use m_constants, only : k_b,pi
   use m_energyfield, only : write_energy_field
   use m_solver_order
   use m_solver_commun
   use m_createspinfile
   use m_update_time
   use m_write_config
   use m_precision
   use m_forces
   use m_velocities
   use m_cell
   use m_io_files_utils
   use m_topo_sd, only: get_charge_map
   use m_topo_commons, only: get_charge, neighbor_Q
   use m_user_info, only: user_info
   use m_FTeff_MD
   use m_vector

   !!!!!!!!!!!!!!!!!!!!!!!
   ! arguments
   !!!!!!!!!!!!!!!!!!!!!!!

   type(lattice), intent(inout)             :: my_lattice
   type(io_parameter), intent(in)           :: io_simu
   type(simulation_parameters), intent(in)  :: ext_param
   type(hamiltonian),intent(inout)          :: H
   ! lattices that are used during the calculations
   type(lattice)    :: lat_1,lat_2

   !!!!!!!!!!!!!!!!!!!!!!!
   ! internal variables
   !!!!!!!!!!!!!!!!!!!!!!!
   real(8) :: E_potential,E_kinetic,E_total
   real(8),allocatable :: masses_motif(:)
   logical :: used(number_different_order_parameters)
   real(8) :: damp_S,damp_F,damp_NVT       ! solid and fluid damping, currently only damp_F is used
   real(8),allocatable :: rand1(:),rand2(:) !tables of random numbers for mkl

   ! arrays to take into account the dynamics
   real(8),allocatable,target       :: V_1(:,:),acc_1(:,:),acc_2(:,:)
   real(8),allocatable,target       :: Feff(:),masses(:),sigma_u(:),sigma_v(:),c_uv(:), delta_vg(:), delta_ug(:)
   real(8),pointer,contiguous       :: Feff_v(:,:),Feff_3(:,:),masses_3(:,:),sigma_u3(:,:),sigma_v3(:,:),c_uv3(:,:), delta_vg3(:,:), delta_ug3(:,:)
   
   ! conversion factor to go from uam.nm/fs^2 to eV/nm
   !real(8), parameter :: uamnmfs_to_eVnm = 9.648526549495E-5  ! to convert the force into an acceleration times a mass
   ! conversion factor to go from uam.nm^2/fs^2 to eV
   !real(8), parameter :: uamnmfs_to_eV = 1.0364277E4   ! to convert the kinetic energy in eV
   
	!conversion of sqrt(kBT/m) in sqrt(eV/uam) to velocity dimension nm/fs
   real(8), parameter ::  convf_v = 9.82269475025328e-03
	! conversion of acceleration in eV/(uma.nm) to nm/fs^2
	real(8), parameter :: convf_a= 96.4853321566533e-006

   ! IOs
   integer :: io_results

   ! dummys
   integer :: N_cell,duration,N_loop,Efreq,gra_freq,j,tag,i,ensemble
   real(8) :: timestep_int,h_int(3),E_int(3),dt
   real(8) :: real_time,Eold,security,Einitial
   real(8) :: kt,ktini,ktfin,Pdy(3),temperature
   logical :: said_it_once,gra_log,io_stochafield,gra_topo
   integer :: dim_mode,iomp !dim_mode of the iterated order parameter
   character(len=100) :: file
   real(8) :: dumy(5),q_plus,q_moins,vortex(3)
   integer,allocatable ::  Q_neigh(:,:)
   real(8) :: time = 0.0d0
   real(8) :: ldc(3)  !Langevin dynamics coefficients
class(ranbase), allocatable :: thermal_noise

   ! prepare the matrices for integration

   call rw_dyna_MD(timestep_int,Efreq,duration,file,damp_S,damp_F,ensemble)
   N_cell=my_lattice%Ncell

   Call my_lattice%used_order(used)
   dim_mode=my_lattice%u%dim_mode
   if(modulo(dim_mode,3)/=0) STOP "dynamics will only work if the considered order-parameter can be described as 3-vector"

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!! Select the propagators and the integrators
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   call select_propagator(ext_param%ktini,N_loop)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!! Create copies of lattice with order-parameter for intermediary states
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   Call my_lattice%copy(lat_1)
   Call my_lattice%copy(lat_2)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!! allocate the element of integrations and associate the pointers to them
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   allocate(Feff(dim_mode*N_cell),source=0.0d0)
   Feff_v(1:dim_mode,1:N_cell)=>Feff
   Feff_3(1:3,1:N_cell*(dim_mode/3))=>Feff

   allocate(V_1(my_lattice%u%dim_mode,N_cell),source=0.0d0)
   allocate(acc_1(my_lattice%u%dim_mode,N_cell),acc_2(my_lattice%u%dim_mode,N_cell),source=0.0d0)
	allocate(rand1(dim_mode*N_cell),rand2(dim_mode*N_cell),source=0.0d0)

   ! get the lattice of the masses
   allocate(masses(dim_mode*N_cell),source=0.0d0)
   Call my_lattice%cell%get_M_phonon(masses_motif)
   do i=1,N_cell
     do j=1,dim_mode
       masses(j+(i-1)*dim_mode)=masses_motif((j-1)/3+1)
     enddo
   enddo
   masses_3(1:3,1:N_cell*(dim_mode/3))=>masses

	allocate(sigma_u(my_lattice%u%dim_mode*N_cell),sigma_v(my_lattice%u%dim_mode*N_cell),c_uv(my_lattice%u%dim_mode*N_cell), &
   			delta_ug(my_lattice%u%dim_mode*N_cell),delta_vg(my_lattice%u%dim_mode*N_cell),source=0.0d0)

   sigma_u3(1:3,1:N_cell*(dim_mode/3))=>sigma_u
   sigma_v3(1:3,1:N_cell*(dim_mode/3))=>sigma_v
   c_uv3(1:3,1:N_cell*(dim_mode/3))=>c_uv
   delta_ug3(1:3,1:N_cell*(dim_mode/3))=>delta_ug
   delta_vg3(1:3,1:N_cell*(dim_mode/3))=>delta_vg

   call initialize_velocities(V_1)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!! start the simulation
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   gra_log=io_simu%io_Xstruct
   io_stochafield=io_simu%io_Tfield
   gra_freq=io_simu%io_frequency
   gra_topo=io_simu%io_topo
   ktini=ext_param%ktini !kBT?
   ktfin=ext_param%ktfin
   kt=ktini
   Eold=100.0d0
   real_time=0.0d0
   Einitial=0.0d0
   h_int=ext_param%H_ext
   E_int=ext_param%E_ext
   said_it_once=.False.
   security=0.0d0
   damp_NVT=0.0d0
   
   !initalize the Langevin parameters [Paterlini and Ferguson, Chem. Phys. 236, 243 (1998)]
   dt=timestep_int !if dt varies the following lines must be moved into the iterations
   ldc(1)=exp(-damp_F*dt) !damp_F is in fs^-1 and damp_F*dt << 1
   ldc(2)=1-(damp_F*dt)/2.0d0+(damp_F*dt)**2/6.0d0-(damp_F*dt)**3/24.0d0+(damp_F*dt)**4/120.0d0     
   ldc(3)= 0.5d0-(damp_F*dt)/6.0d0+(damp_F*dt)**2/24.0d0-(damp_F*dt)**3/120.0d0  

	if (kt.gt.1.0d-8) then !warning: these depend on dt
		if(damp_F.eq.0.0d0) stop 'for T!=0 damp_F should not be 0'
		sigma_v(:)=  convf_v *sqrt( kT/masses(:) * ( 1.0d0-exp(-2.0d0*damp_F*dt) ) )  	!bivariate gaussian distrib variance 1,  in unit of the velocities: nm/fs
  		sigma_u(:)=  convf_v *sqrt( dt*kT/masses(:)/(damp_F) * ( 2.0d0 - (3.0d0 - 4.0d0*exp(-damp_F*dt) + exp(-2.0d0*damp_F*dt)) / (damp_F*dt) ) )  !bivariate gaussian distrib variance 2, in unit of u: nm
   	c_uv(:)= convf_v**2 * kT/masses(:) *1.0d0/(sigma_u(:)*sigma_v(:)*damp_F) * (1.0d0 - exp(-damp_F*dt) )**2  !correlation coefficient, dimensionless
	endif

	write(6,'(a,3(2x,E20.12E3))') 'Warning: thermostat implementation in MD for constant timestep only, change it for variable timestep.'

	
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! initialize the simulation
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Call my_lattice%copy_val_to(lat_1)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!! Prepare Q_neigh for topological neighbor calculation
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call user_info(output_unit,time,'topological operators',.false.)
    Call neighbor_Q(my_lattice,Q_neigh)
    call user_info(output_unit,time,'done',.true.)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! initialize the random number in case of thermal noise
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (kt.gt.1.0d-8) then
       call get_ran_type(thermal_noise)
       call thermal_noise%init_base(dim_mode*N_cell)
    endif

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! initialize Energy
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!

    E_potential=H%energy(my_lattice)/real(N_cell,8)
    E_kinetic=0.5d0*sum(masses_3*V_1**2)/real(N_cell,8)/convf_v**2 !convert to eV
    Eold=E_potential+E_kinetic

    write(output_unit,'(a,3(2x,E20.12E3))') 'Initial potential, kinetic and Total Energy (eV)',E_potential,E_kinetic,Eold

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! do we use the update timestep
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call init_update_time('input')

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! prepare the output
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!

    io_results=open_file_write('EM.dat')
    Write(io_results,'(23(a,2x))') '# 1:step','2:real_time','3:E_av','4:U', &
         &  '5:Ux','6:Uy','7:Uz','8:vorticity','9:vx', &
         &  '10:vy','11:vz','12:qeuler','13:q+','14:q-','15:T=', &
         &  '16:Epot=','17:Ek=','18:Ex','19:Ey=','20:Ez=','21:Hx','22:Hy=','23:Hz='


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! initialize acceleration
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!get forces on the phonon lattice
	Call H%get_eff_field(lat_1,Feff,5) !Feff here is eV/nm
	acc_1=convf_a*Feff_3/masses_3 !in nm/fs^2 a(t)
                    
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! beginning of the simulation
	do j=1,duration
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		call truncate(lat_1,used)
      !dt=timestep_int

      !!!!!!!!!!!!!!!!!!!
      ! Stochastic Velocity Verlet algorithm
		!!!!!!!!!!!!!!!!!!!
   	
      !if non zero temperature: draw random numbers delta_vg, delta_ug           
      if (kt.gt.1.0d-8) call draw_stocha_integrals(thermal_noise,rand1,rand2,sigma_u,sigma_v,c_uv,delta_ug,delta_vg)
 		   
		!update positions
       lat_2%u%modes_3= lat_1%u%modes_3  + ldc(2)*V_1*dt + ldc(3)*acc_1*dt**2 + delta_ug3     !u(t+dt) in nm
      
      !update accelerations
      Call H%get_eff_field(lat_2,Feff,5)
      acc_2=convf_a*Feff_3/masses_3 !a(t+dt)

		!update velocities
       V_1=ldc(1)*V_1 + (ldc(2)-ldc(3))*acc_1*dt + ldc(3)*acc_2*dt + delta_vg3   !v(t+dt) 
       
       !swap accelerations
       acc_1=acc_2

      !!!!!!!!!!!!!!!!!!!
      ! end Verlet algorithm
      !!!!!!!!!!!!!!!!!!!

      Call lat_2%u%copy_val(lat_1%u)
      call truncate(lat_1,used)
		
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!! plotting with graphical frequency
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(mod(j-1,gra_freq)==0)then
            tag=j/gra_freq

            E_potential=H%energy(lat_1)/real(N_cell,8)
            E_kinetic=0.5d0*sum(masses_3*V_1**2)/real(N_cell,8)/convf_v**2 !convert to eV
            E_total=E_potential+E_kinetic
            temperature=2.0d0*E_kinetic/3.0d0/k_b

            Pdy=sum(lat_1%u%modes_3,2)/real(N_cell,8) !sums of the displacements over all atoms in unit cell

            dumy=get_charge(lat_1,Q_neigh)
            q_plus=dumy(1)/pi/4.0d0
            q_moins=dumy(2)/pi/4.0d0
            vortex=dumy(3:5)/3.0d0/sqrt(3.0d0)

            if(gra_log) then
                call CreateSpinFile(tag,lat_1%u)
                Call write_config(tag,lat_1)

!                call write_netcdf('test',lat_1,real_time)
                write(output_unit,'(a,3x,I10)') 'wrote phonon configuration and povray file number',tag
                write(output_unit,'(a,3x,f14.6,3x,a,3x,I10)') 'real time in ps',real_time/1000.0d0,'iteration',j
            endif

            Write(io_results,'(I6,21(E20.12E3,2x),E20.12E3)') j,real_time,E_total, &
             &   norm2(Pdy),Pdy,norm2(vortex),vortex,q_plus+q_moins,q_plus,q_moins, &
             &   temperature,E_potential,E_kinetic,E_int,H_int

            if (io_simu%io_Force)& !call forces(tag,lat_1%ordpar%all_l_modes,my_lattice%dim_mode,my_lattice%areal)
                & ERROR STOP "FORCES HAVE TO BE REIMPLEMENTED"

            write(output_unit,'(a,4(2x,E20.12E3))') 'E_pot, E_k, E_tot (eV) and T (K)',E_potential,E_kinetic,E_total,temperature
            write(output_unit,'(a,2x,I8,2x,E20.12E3,/)') 'step and time (in fs)',j,real_time+timestep_int
        endif

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! update timestep
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!

!        call update_time(timestep_int,Feff_v,damping)tail 

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!! end of a timestep
        real_time=real_time+timestep_int !increment time counter

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   enddo
   ! end of the simulation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call close_file('EM.dat',io_results)

Call lat_1%copy_val_to(my_lattice)

end subroutine

end module
