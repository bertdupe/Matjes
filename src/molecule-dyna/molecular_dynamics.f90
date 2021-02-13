!
!subroutine that carries out molecular dynamics

module m_molecular_dynamics

implicit none
contains

subroutine molecular_dynamics(my_lattice,motif,io_simu,ext_param,Hams)
   use m_derived_types, only : lattice,io_parameter,number_different_order_parameters,simulation_parameters
   use m_energy_output_contribution, only:Eout_contrib_init, Eout_contrib_write
   use m_constants, only : k_b
   use m_energyfield, only : write_energy_field
   use m_solver_order
   use m_H_public
   use m_solver_commun
   use m_eff_field, only: get_eff_field
   use m_createspinfile
   use m_update_time
   use m_write_config
   use m_precision
   use m_forces
   use m_velocities
   use m_cell

   !!!!!!!!!!!!!!!!!!!!!!!
   ! arguments
   !!!!!!!!!!!!!!!!!!!!!!!

   type(lattice), intent(inout) :: my_lattice
   type(t_cell), intent(in) :: motif
   type(io_parameter), intent(in) :: io_simu
   type(simulation_parameters), intent(in) :: ext_param
   class(t_H), intent(in) :: Hams(:)
   ! lattices that are used during the calculations
   type(lattice)                         :: lat_1,lat_2

   !!!!!!!!!!!!!!!!!!!!!!!
   ! internal variables
   !!!!!!!!!!!!!!!!!!!!!!!
   real(8) :: E_potential,E_kinetic,E_total
   real(8),allocatable :: masses_motif(:)
   logical :: used(number_different_order_parameters)

   ! arrays to take into account the dynamics
   real(8),allocatable,target       :: Du(:,:,:),Du_int(:,:),V_1(:,:),V_2(:,:),acceleration(:,:),masses(:,:)
   real(8),allocatable,target       :: Feff(:)
   real(8),pointer,contiguous       :: Feff_v(:,:),Feff_3(:,:)
   real(8),pointer,contiguous       :: Du_3(:,:,:),Du_int_3(:,:)

   ! dummys
   integer :: N_cell,duration,N_loop,Efreq,gra_freq,j,i_loop,tag,i
   real(8) :: timestep_int,h_int(3),E_int(3),dt
   real(8) :: real_time,Eold,security,Einitial
   real(8) :: kt,ktini,ktfin,Pdy(3)
   logical :: said_it_once,gra_log,io_stochafield,gra_topo
   integer :: dim_mode !dim_mode of the iterated order parameter
   character(len=100) :: file

   ! prepare the matrices for integration
   call rw_dyna_MD(timestep_int,Efreq,duration,file)
   N_cell=product(my_lattice%dim_lat)
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

   E_potential=energy_all(Hams,my_lattice)

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
   allocate(masses(my_lattice%u%dim_mode,N_cell),source=0.0d0)
   Call motif%get_M_phonon(masses_motif)
   do i=1,N_cell
     do j=1,my_lattice%u%dim_mode
       masses(j,i)=masses_motif((j-1)/3+1)
     enddo
   enddo

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

    E_potential=energy_all(Hams,my_lattice)/real(N_cell,8)
    E_kinetic=0.5d0*sum(masses*V_1**2)/real(N_cell,8)
    Eold=E_potential+E_kinetic

    write(6,'(a,3(2x,E20.12E3))') 'Initial potential, kinetic and Total Energy (eV)',E_potential,E_kinetic,Eold

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! do we use the update timestep
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call init_update_time('input')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! beginning of the simulation
    do j=1,duration
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call truncate(lat_1,used)
        dt=timestep_int

        !
        ! version 1 of the integration
        !

        !get forces on the phonon lattice
        Call get_eff_field(Hams,lat_1,Feff,5)

        acceleration=Feff_3/masses

        V_2=acceleration*dt+V_1

        lat_2%u%modes_3=V_2*dt+lat_1%u%modes_3











!        !
!        ! loop over the integration order
!        !
!        do i_loop=1,N_loop
!          !get actual dt from butchers table
!          dt=get_dt_mode(timestep_int,i_loop)
!
!        !update phonon
!          !get effective field on magnetic lattice
!          Call get_eff_field(Hams,lat_1,Feff,5)
!          !do integration
!          Du_3(:,:,i_loop)=Feff_3/1.0d0
!          pause
!          Call get_Dmode_int(Du,i_loop,N_loop,Du_int)
!          lat_2%u%modes_3=get_integrator_field(my_lattice%u%modes_3,Du_int_3,dt)
!        !copy phonon
!          Call lat_2%u%copy_val(lat_1%u)
!        enddo
        !!!!!!!!!!!!!!! copy the final configuration in my_lattice
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        V_1=V_2
        Call lat_2%u%copy_val(my_lattice%u)
        call truncate(my_lattice,used)

        E_potential=energy_all(Hams,my_lattice)/real(N_cell,8)
        E_kinetic=0.5d0*sum(masses*V_2**2)/real(N_cell,8)
        E_total=E_potential+E_kinetic

        !if (dabs(check(2)).gt.1.0d-8) call get_temp(security,check,kt)

        if (mod(j-1,Efreq).eq.0) then
            Pdy=sum(my_lattice%u%modes_3,2)/real(N_cell,8) !sums over all atoms in unit cell
        endif

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!! plotting with graphical frequency
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(mod(j-1,gra_freq)==0)then
            tag=j/gra_freq

            if(gra_log) then
                call CreateSpinFile(tag,my_lattice%u)
                Call write_config(tag,my_lattice)
                write(6,'(a,3x,I10)') 'wrote phonon configuration and povray file number',tag
                write(6,'(a,3x,f14.6,3x,a,3x,I10)') 'real time in ps',real_time/1000.0d0,'iteration',j
            endif

            if (io_simu%io_Force)& !call forces(tag,lat_1%ordpar%all_l_modes,my_lattice%dim_mode,my_lattice%areal)
                & ERROR STOP "FORCES HAVE TO BE REIMPLEMENTED"

        endif

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! update timestep
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!

!        call update_time(timestep_int,Feff_v,damping)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!! end of a timestep
        real_time=real_time+timestep_int !increment time counter

        write(6,'(a,3(2x,E20.12E3))') 'Potential, kinetic and Total Energy (eV)',E_potential,E_kinetic,E_total
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   enddo
   ! end of the simulation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine

end module
