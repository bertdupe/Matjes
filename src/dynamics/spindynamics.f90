module m_spindynamics
use m_derived_types, only : lattice,number_different_order_parameters
use m_derived_types, only : io_parameter,simulation_parameters
use m_hamiltonian_collection, only: hamiltonian
use m_dyna_io, only: dyna_input
use mpi_basic, only: mpi_type

implicit none
contains

subroutine spindynamics(lat,io_simu,ext_param,Ham,Ham_res,comm)
    !wrapper to first initialize all spin-dynamics parameters on the necessary threads
    use m_H_public,only: t_H
    use m_dipolar_fft, only: get_dip
    type(lattice), intent(inout)                :: lat
    type(io_parameter), intent(inout)           :: io_simu
    type(simulation_parameters), intent(inout)  :: ext_param
    class(t_H),intent(in),allocatable           :: Ham(:),Ham_res(:)
    type(mpi_type),intent(in)                   :: comm

    type(hamiltonian)       :: H_res,H
    type(dyna_input)        :: io_dyn

    if(comm%ismas)then
        Call H%init_H_cp(Ham)   !later change to move, as certain the result is the same (do it even in setup_simu?)
        if(allocated(Ham_res))  Call H_res%init_H_cp(Ham_res)
        Call get_dip(H%dip,lat)
        H%NH_total=H%NH_total+1
        allocate(H%dip_H(3*lat%nmag,lat%ncell))
        if(allocated(Ham_res))then
            Call get_dip(H_res%dip,lat)
            H_res%NH_total=H_res%NH_total+1
            allocate(H_res%dip_H(3*lat%nmag,lat%ncell))
        endif
    endif

    !fill Hamiltonian which allows the distribution
    Call H%distribute(comm)
    Call H_res%distribute(comm)

    Call lat%bcast(comm)
    Call io_simu%bcast(comm)
    Call ext_param%bcast(comm)

    Call io_dyn%read_file(comm=comm)

    Call spindynamics_run(lat,io_dyn,io_simu,ext_param,H,H_res,comm)
end subroutine

subroutine spindynamics_run(mag_lattice,io_dyn,io_simu,ext_param,H,H_res,comm)
    use m_topo_commons, only: get_charge, neighbor_Q
    use m_update_time, only: update_time,init_update_time
    use m_constants, only : pi,k_b
    use m_energyfield, only : write_energy_field 
    use m_createspinfile, only: CreateSpinFile
    use m_user_info, only: user_info
    use m_excitations, only: update_ext_EM_fields, update_EMT_of_r,set_excitations
    use m_solver_commun, only: get_integrator_field, get_propagator_field,select_propagator
    use m_topo_sd, only: get_charge_map
    use m_solver_order,only: get_dt_mode
    use m_tracker, only: init_tracking,plot_tracking
    use m_print_Beff, only: print_Beff
    use m_precision, only: truncate
    use m_H_public, only: t_H, energy_all
    use m_eff_field, only :get_eff_field
    use m_write_config, only: write_config
    use m_energy_output_contribution, only:Eout_contrib_init, Eout_contrib_write
    use m_solver_order,only : get_Dmode_int
    use m_dyna_io, only: rw_dyna
!$  use omp_lib
    
    ! input
    type(lattice), intent(inout)    :: mag_lattice
    type(dyna_input),intent(in)     :: io_dyn
    type(io_parameter), intent(in)  :: io_simu
    type(simulation_parameters), intent(in) :: ext_param
    type(hamiltonian),intent(inout)     ::  H,H_res
    type(mpi_type),intent(in)       :: comm
    ! internal
    logical :: gra_log,io_stochafield
    integer :: i,j,i_loop,input_excitations
    ! lattices that are used during the calculations
    type(lattice)                         :: lat_1,lat_2
    
    !intermediate values for dynamics
    real(8),allocatable,target              :: Dmag(:,:,:),Dmag_int(:,:)
    real(8),allocatable,dimension(:),target :: Beff(:)
    real(8),pointer,contiguous              :: Beff_v(:,:)
    real(8),pointer,contiguous              :: Beff_3(:,:)
    real(8),pointer,contiguous              :: Dmag_3(:,:,:),Dmag_int_3(:,:)
    
    ! dummys
    real(8) :: q_plus,q_moins,vortex(3),Mdy(3),Edy,Eold,dt
    real(8) :: check(2),test_torque,Einitial,ave_torque
    real(8) :: dumy(5),security(2)
    real(8) :: timestep_int,real_time,h_int(3),E_int(3)
    real(8) :: kt,ktini,ktfin
    real(8) :: time
    integer :: iomp, N_cell, N_loop, tag
    !integer :: io_test
    !! switch that controls the presence of magnetism, electric fields and magnetic fields
    logical :: i_excitation
    logical :: used(number_different_order_parameters)
    ! dumy
    logical :: said_it_once,gra_topo
    integer,allocatable ::  Q_neigh(:,:)
    integer :: io_Eout_contrib
    integer :: dim_mode !dim_mode of the iterated order parameter

    real(8),allocatable ::  Beff_tmp(:)


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!! Initialize all data necessary on all threads
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    dim_mode=mag_lattice%M%dim_mode
    N_cell=mag_lattice%Ncell
    Call mag_lattice%used_order(used)

    ! Create copies of lattice with order-parameter for intermediary states
    Call mag_lattice%copy(lat_1) 
    allocate(Beff(dim_mode*N_cell),source=0.0d0)
    allocate(Beff_tmp(dim_mode*N_cell),source=0.0d0)

    !generalize the following
    call select_propagator(ext_param%ktini,N_loop)

    if(comm%ismas)then
        time=0.0d0
        input_excitations=0
        if(modulo(dim_mode,3)/=0) STOP "dynamics will only work if the considered order-parameter can be described as 3-vector"
        
        OPEN(7,FILE='EM.dat',action='write',status='replace',form='formatted')
              Write(7,'(20(a,2x))') '# 1:step','2:real_time','3:E_av','4:M', &
             &  '5:Mx','6:My','7:Mz','8:vorticity','9:vx', &
             &  '10:vy','11:vz','12:qeuler','13:q+','14:q-','15:T=', &
             &  '16:Tfin=','17:Ek=','18:Hx','19:Hy=','20:Hz='
        
        ! check the convergence
        open(8,FILE='convergence.dat',action='write',form='formatted')

        if(io_simu%io_energy_cont)then
            if(io_simu%io_energy_detail)then
                Call Eout_contrib_init(H_res,io_Eout_contrib)
            else
                Call Eout_contrib_init(H,io_Eout_contrib)
            endif
        endif
        
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!! Prepare Q_neigh for topological neighbor calculation
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        call user_info(6,time,'topological operators',.false.)
        Call neighbor_Q(mag_lattice,Q_neigh)
        call user_info(6,time,'done',.true.)
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!! allocate the element of integrations and associate the pointers to them
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        Call mag_lattice%copy(lat_2) 

        allocate(Dmag(mag_lattice%M%dim_mode,N_cell,N_loop),source=0.0d0) 
        allocate(Dmag_int(mag_lattice%M%dim_mode,N_cell),source=0.0d0) 
        Dmag_3(1:3,1:N_cell*(dim_mode/3),1:N_loop)=>Dmag
        Dmag_int_3(1:3,1:N_cell*(dim_mode/3))=>Dmag_int

        Beff_v(1:dim_mode,1:N_cell)=>Beff
        Beff_3(1:3,1:N_cell*(dim_mode/3))=>Beff
        
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!! start the simulation
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        gra_log=io_simu%io_Xstruct
        io_stochafield=io_simu%io_Tfield
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
        timestep_int=io_dyn%timestep

        Call mag_lattice%copy_val_to(lat_1) !necessary due to loop structure
    endif
    
    Edy=H%energy(mag_lattice)

    if(comm%ismas)then
        write(6,'(a,2x,E20.12E3)') 'Initial Total Energy (eV)',Edy/real(N_cell,8)
    
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! part of the excitations
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call set_excitations('input',i_excitation,input_excitations)
    
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! do we use the update timestep
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call init_update_time('input')
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! initialize the different torques
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !call get_torques('input')
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! check if a magnetic texture should be tracked
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (io_simu%io_tracker) call init_tracking(mag_lattice)
    endif
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! beginning of the simulation
    do j=1,io_dyn%duration
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !   call init_temp_measure(check,check1,check2,check3)
        if(comm%ismas)then
            call truncate(lat_1,used)
            Edy=0.0d0
            ave_torque=0.0d0
            test_torque=0.0d0
            dt=timestep_int
            !why is this outside of the integration order loop? time changes there
            call update_ext_EM_fields(real_time,check)
        endif
       
       !
       ! loop over the integration order
       !
       do i_loop=1,N_loop
           !get actual dt from butchers table
            if(comm%ismas)then
                dt=get_dt_mode(timestep_int,i_loop)
                
                ! loop that get all the fields
                if (i_excitation) then
                    do iomp=1,N_cell
                        !smarter to just copy relevant order parameters around, or even point all to the same array
                        call update_EMT_of_r(iomp,mag_lattice)
                        call update_EMT_of_r(iomp,lat_1)
                    enddo
                endif
            endif

           !get effective field on magnetic lattice
            if(comm%Np>1) Call lat_1%bcast_val(comm)
           Call H%get_eff_field(lat_1,Beff,1,Beff_tmp)

            if(comm%ismas)then
                !do integration
                ! Be carefull the sqrt(dt) is not included in BT_mag(iomp),D_T_mag(iomp) at this point. It is included only during the integration
                Call get_propagator_field(Beff_3,io_dyn%damping,lat_1%M%modes_3,Dmag_3(:,:,i_loop))
                Call get_Dmode_int(Dmag,i_loop,N_loop,Dmag_int)
                lat_2%M%modes_3=get_integrator_field(mag_lattice%M%modes_3,Dmag_int_3,dt)
                !copy magnetic texture to 1 
                Call lat_2%M%copy_val(lat_1%M)
            endif
        enddo
        !!!!!!!!!!!!!!! copy the final configuration in my_lattice
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(comm%ismas)then
            Call lat_2%M%copy_val(mag_lattice%M)
            call truncate(mag_lattice,used)
        endif
       
        if(comm%Np>1) Call mag_lattice%bcast_val(comm)
        Edy=H%energy(mag_lattice)/real(N_cell,8)
       
        if(comm%ismas)then
            if (mod(j-1,io_dyn%Efreq).eq.0) then
                !get values to plot (Mavg,topo)
                Mdy=sum(mag_lattice%M%modes_3,2)/real(N_cell,8)  !sums over all magnetic atoms in unit cell ( not sure if this is wanted)
                if(mag_lattice%nmag>1) ERROR STOP "TOPO CHARGE WILL PROBABLY NOT WORK FOR nmag>1"
                dumy=get_charge(lat_1,Q_neigh)
                q_plus=dumy(1)/pi/4.0d0
                q_moins=dumy(2)/pi/4.0d0
                vortex=dumy(3:5)/3.0d0/sqrt(3.0d0)
                !write data files
                Write(7,'(I12,18(E20.12E3,2x),E20.12E3)') j,real_time,Edy, &
                 &   norm2(Mdy),Mdy,norm2(vortex),vortex,q_plus+q_moins,q_plus,q_moins, &
                 &   kT/k_B,(security(i),i=1,2),H_int
                write(8,'(I10,3x,3(E20.12E3,3x))') j,Edy,test_torque,ave_torque
                if(io_simu%io_energy_cont)then
                    if(io_simu%io_energy_detail)then
                        Call Eout_contrib_write(H_res,j,real_time,mag_lattice,io_Eout_contrib)
                    else
                        Call Eout_contrib_write(H,j,real_time,mag_lattice,io_Eout_contrib)
                    endif

                endif
            endif
            
            ! security in case of energy increase in SD and check for convergence
            if (((io_dyn%damping*(Edy-Eold).gt.1.0d-10).or.(io_dyn%damping*(Edy-Einitial).gt.1.0d-10)).and.(kt.lt.1.0d-10).and.(.not.said_it_once)) then
                write(6,'(a)') 'WARNING: the total energy or torque is increasing for non zero damping'
                write(6,'(a)') 'this is not allowed by theory'
                write(6,'(a)') 'please reduce the time step'
                said_it_once=.True.
            endif
            Eold=Edy
        endif
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!! plotting with graphical frequency
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(mod(j-1,io_simu%io_frequency)==0)then
            tag=j/io_simu%io_frequency

            if (io_simu%io_Energy_Distrib)then
                if(comm%Np>1) Call mag_lattice%bcast_val(comm)
                if(io_simu%io_Energy_detail)then
                    Call write_energy_field(tag,H_res,mag_lattice,1,comm%ismas) !1 for magnetic field
                else
                    Call write_energy_field(tag,H    ,mag_lattice,1,comm%ismas) !1 for magnetic field
                endif
            endif

            if(comm%ismas)then
                if (io_simu%io_Beff) call print_Beff(tag,Beff_v)

                if (io_simu%io_tracker)then
                    ERROR STOP "plot_tracking is not implemented with the new lattce and Hamiltonian"
                    !call plot_tracking(j/io_simu%io_frequency,lat_1,Hams)
                endif

                if(gra_log) then
                    call CreateSpinFile(tag,mag_lattice%M)
                    Call write_config(tag,mag_lattice) 
                    write(6,'(a,3x,I10)') 'wrote Spin configuration and povray file number',tag
                    write(6,'(a,3x,f14.6,3x,a,3x,I10)') 'real time in ps',real_time/1000.0d0,'iteration',j
                endif
                if(gra_topo) Call get_charge_map(tag,mag_lattice,Q_neigh)
       
                if (io_simu%io_Force) &!call forces(tag,lat_1%ordpar%all_l_modes,mag_lattice%dim_mode,mag_lattice%areal)
                    & ERROR STOP "FORCES HAVE TO BE REIMPLEMENTED"
       
                if(io_simu%io_fft_Xstruct) &!call plot_fft(mag_lattice%ordpar%all_l_modes,-1.0d0,mag_lattice%areal,mag_lattice%dim_lat,mag_lattice%periodic,mag_lattice%dim_mode,tag)
                    & ERROR STOP "PLOT FFT HAS TO BE REIMPLEMENTED"
            endif
        endif
             
        if(comm%ismas)then
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! update timestep
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
            call update_time(timestep_int,Beff_v,io_dyn%damping)
            
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! reinitialize T variables
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!
            check(1)=0.0d0
            check(2)=0.0d0
            
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!! end of a timestep
            real_time=real_time+timestep_int !increment time counter
        endif
    enddo 
    if(comm%ismas)then 
        !!!!!!!!!!!!!!! end of iteration
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        close(7)
        close(8)
        if(io_simu%io_energy_cont) close(io_Eout_contrib)

        write(6,'(a,2x,E20.12E3)') 'Final Total Energy (eV)',Edy
        if ((dabs(check(2)).gt.1.0d-8).and.(kt/k_B.gt.1.0d-5)) then
            write(6,'(a,2x,f16.6)') 'Final Temp (K)', check(1)/check(2)/2.0d0/k_B
            write(6,'(a,2x,f14.7)') 'Kinetic energy (meV)', (check(1)/check(2)/2.0d0-kT)/k_B*1000.0d0
        else
            write(6,'(a)') 'the temperature measurement is not possible'
        endif
    endif

end subroutine 

end module
