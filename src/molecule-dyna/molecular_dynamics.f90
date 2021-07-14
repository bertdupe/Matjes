!
!subroutine that carries out molecular dynamics
module m_molecular_dynamics
use m_derived_types, only : lattice,number_different_order_parameters
use m_derived_types, only : io_parameter,simulation_parameters
use m_hamiltonian_collection, only: hamiltonian
use mpi_basic, only: mpi_type
use m_H_public
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
          write(error_unit,'(A)')   "         Only one thead will to anything."
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
   real(8) :: damp_S,damp_F       ! solid and fluid damping

   ! arrays to take into account the dynamics
   real(8),allocatable,target       :: Du(:,:,:),Du_int(:,:),V_1(:,:),V_2(:,:),acceleration(:,:)
   real(8),allocatable,target       :: Feff(:),FT_eff(:),masses(:)
   real(8),pointer,contiguous       :: Feff_v(:,:),Feff_3(:,:),FT_eff_3(:,:),masses_3(:,:)
   real(8),pointer,contiguous       :: Du_3(:,:,:),Du_int_3(:,:)

   ! conversion factor to go from uam.nm/fs^2 to eV/nm
   real(8), parameter :: uamnmfs_to_eVnm = 9.648526549495E-5  ! to convert the force into an accelaration times a mass
   ! conversion factor to go from uam.nm^2/fs^2 to eV
   real(8), parameter :: uamnmfs_to_eV = 1.0364277E4   ! to convert the kinetic energy in eV

   ! IOs
   integer :: io_results

   ! dummys
   integer :: N_cell,duration,N_loop,Efreq,gra_freq,j,tag,i,ensemble
   real(8) :: timestep_int,h_int(3),E_int(3),dt
   real(8) :: real_time,Eold,security,Einitial
   real(8) :: kt,ktini,ktfin,Pdy(3),temperature
   logical :: said_it_once,gra_log,io_stochafield,gra_topo
   integer :: dim_mode !dim_mode of the iterated order parameter
   character(len=100) :: file
   real(8) :: dumy(5),q_plus,q_moins,vortex(3)
   integer,allocatable ::  Q_neigh(:,:)
   real(8) :: time = 0.0d0

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

   allocate(Feff(my_lattice%u%dim_mode*N_cell),source=0.0d0)
   Feff_v(1:my_lattice%u%dim_mode,1:N_cell)=>Feff
   Feff_3(1:3,1:N_cell*(dim_mode/3))=>Feff

   allocate(Du(my_lattice%u%dim_mode,N_cell,N_loop),source=0.0d0)
   allocate(Du_int(my_lattice%u%dim_mode,N_cell),source=0.0d0)
   Du_3(1:3,1:N_cell*(dim_mode/3),1:N_loop)=>Du
   Du_int_3(1:3,1:N_cell*(dim_mode/3))=>Du_int

   allocate(V_1(my_lattice%u%dim_mode,N_cell),source=0.0d0)
   allocate(V_2(my_lattice%u%dim_mode,N_cell),source=0.0d0)
   allocate(acceleration(my_lattice%u%dim_mode,N_cell),source=0.0d0)

   ! get the lattice of the masses
   allocate(masses(my_lattice%u%dim_mode*N_cell),source=0.0d0)
   Call my_lattice%cell%get_M_phonon(masses_motif)
   do i=1,N_cell
     do j=1,my_lattice%u%dim_mode
       masses(j+(i-1)*my_lattice%u%dim_mode)=masses_motif((j-1)/3+1)
     enddo
   enddo
   masses_3(1:3,1:N_cell*(dim_mode/3))=>masses

   allocate(FT_eff(my_lattice%u%dim_mode*N_cell),source=0.0d0)
   FT_eff_3(1:3,1:N_cell*(dim_mode/3))=>FT_eff

   call initialize_velocities(V_1)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!! start the simulation
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   gra_log=io_simu%io_Xstruct
   io_stochafield=io_simu%io_Tfield
   gra_freq=io_simu%io_frequency
   gra_topo=io_simu%io_topo
   ktini=ext_param%ktini
   ktfin=ext_param%ktfin
   kt=ktini
   Eold=100.0d0
   real_time=0.0d0
   Einitial=0.0d0
   h_int=ext_param%H_ext
   E_int=ext_param%E_ext
   said_it_once=.False.
   security=0.0d0

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! initialize the simulation
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Call my_lattice%copy_val_to(lat_1)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!! Prepare Q_neigh for topological neighbor calculation
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call user_info(6,time,'topological operators',.false.)
    Call neighbor_Q(my_lattice,Q_neigh)
    call user_info(6,time,'done',.true.)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! initialize Energy
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!

    E_potential=H%energy(my_lattice)/real(N_cell,8)
    E_kinetic=uamnmfs_to_eV*0.5d0*sum(masses_3*V_1**2)/real(N_cell,8)
    Eold=E_potential+E_kinetic

    write(6,'(a,3(2x,E20.12E3))') 'Initial potential, kinetic and Total Energy (eV)',E_potential,E_kinetic,Eold

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
    ! beginning of the simulation
    do j=1,duration
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call truncate(lat_1,used)
        dt=timestep_int

        !!!!!!!!!!!!!!!!!!!
        ! Verlet algorithm
        !!!!!!!!!!!!!!!!!!!

           !get forces on the phonon lattice
           Call H%get_eff_field(lat_1,Feff,5)
           if (kt.gt.1.0d-8) call Bolzmann(kt,damp_F,masses,FT_eff)
           acceleration=(uamnmfs_to_eVnm*(Feff_3+FT_eff_3)-damp_F*V_1)/masses_3

           V_2=acceleration*dt/2.0d0+V_1      ! ( v of t+dt/2  )
           lat_2%u%modes_3=V_2*dt+lat_1%u%modes_3       ! ( r of t+dt  )

           !get forces on the phonon lattice

           Call H%get_eff_field(lat_2,Feff,5)
           acceleration=(uamnmfs_to_eVnm*Feff_3-damp_F*V_2)/masses_3
           V_1=acceleration*dt/2.0d0+V_2      ! ( v of t+dt  )

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
            E_kinetic=uamnmfs_to_eV*0.5d0*sum(masses_3*V_1**2)/real(N_cell,8)
            E_total=E_potential+E_kinetic
            temperature=2.0d0*E_kinetic/3.0d0/real(N_cell,8)/k_b

            Pdy=sum(lat_1%u%modes_3,2)/real(N_cell,8) !sums of the displacements over all atoms in unit cell

            dumy=get_charge(lat_1,Q_neigh)
            q_plus=dumy(1)/pi/4.0d0
            q_moins=dumy(2)/pi/4.0d0
            vortex=dumy(3:5)/3.0d0/sqrt(3.0d0)

            if(gra_log) then
                call CreateSpinFile(tag,lat_1%u)
                Call write_config(tag,lat_1)
                write(6,'(a,3x,I10)') 'wrote phonon configuration and povray file number',tag
                write(6,'(a,3x,f14.6,3x,a,3x,I10)') 'real time in ps',real_time/1000.0d0,'iteration',j
            endif

            Write(io_results,'(I6,21(E20.12E3,2x),E20.12E3)') j,real_time,E_total, &
             &   norm2(Pdy),Pdy,norm2(vortex),vortex,q_plus+q_moins,q_plus,q_moins, &
             &   temperature,E_potential,E_kinetic,E_int,H_int

            if (io_simu%io_Force)& !call forces(tag,lat_1%ordpar%all_l_modes,my_lattice%dim_mode,my_lattice%areal)
                & ERROR STOP "FORCES HAVE TO BE REIMPLEMENTED"

            write(6,'(a,4(2x,E20.12E3))') 'E_pot, E_k, E_tot (eV) and T (K)',E_potential,E_kinetic,E_total,temperature
            write(6,'(a,2x,I8,2x,E20.12E3,/)') 'step and time (in fs)',j,real_time+timestep_int
        endif

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! update timestep
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!

!        call update_time(timestep_int,Feff_v,damping)

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
