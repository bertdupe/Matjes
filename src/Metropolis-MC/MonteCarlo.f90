module m_montecarlo
implicit none
contains
!
! Routine that does the Monte Carlo (and not the parallel tempering)
!
!subroutine MonteCarlo(my_lattices,mag_motif,io_simu,gra_topo,ext_param)
subroutine montecarlo(my_lattice,io_simu,ext_param,Hams)
    use m_constants, only : k_b,pi
    use m_vector, only : norm
    use m_derived_types, only : lattice,t_cell,io_parameter,simulation_parameters
    use m_basic_types, only : vec_point
    use m_rw_MC
    use m_topo_commons, only : neighbor_Q,get_charge
    use m_convert
    use m_io_files_utils
    use m_operator_pointer_utils
    use m_eval_Beff
    use m_MCstep
    use m_H_public
    use m_relaxation
    use m_write_config

    use m_Corre
    use m_check_restart
    use m_qorien
    use m_check_restart
    use m_average_MC
    use m_write_spin
    use m_set_temp
    use m_topocharge_all
    use m_createspinfile
    use m_derived_types

    use m_input_types,only: MC_input

    type(lattice), intent(inout) :: my_lattice
    type(io_parameter), intent(in) :: io_simu
    type(simulation_parameters), intent(in) :: ext_param
    class(t_H), intent(in) :: Hams(:)
    
    !!!!!!!!!!!!!!!!!!!!!!!
    ! internal variables
    !!!!!!!!!!!!!!!!!!!!!!!
    ! slope of the MC
    integer :: i_relax,n_kT,n_MC,N_cell,size_table
    ! variable for the temperature
    real(kind=8) :: kT,kTini,kTfin
    ! variables that being followed during the simulation
    real(kind=8) :: qeulerp,qeulerm,vortex(3),magnetization(3),E_total,dumy(5)
    ! contribution of the different energy parts
    real(kind=8) :: E_decompose(8)
    ! thermodynamical quantities
    real(kind=8),allocatable :: C_av(:),chi_M(:,:),chi_Q(:,:)
    ! errors on the different quantities
    real(kind=8),allocatable :: E_err_av(:),M_err_av(:,:)
    ! sums
    real(kind=8),allocatable :: M_sq_sum_av(:,:),E_sum_av(:),E_sq_sum_av(:),Q_sq_sum_av(:),Qp_sq_sum_av(:),Qm_sq_sum_av(:)
    ! energy and so on
    real(kind=8),allocatable :: E_av(:),qeulerp_av(:),qeulerm_av(:),kt_all(:),M_av(:,:)
    real(kind=8),allocatable :: M_sum_av(:,:),vortex_av(:,:),chi_l(:,:)
    ! variable for the convergence of the MC
    real(kind=8) :: acc,rate,tries,cone,nb
    integer :: i
    type(MC_input)  :: io_MC
    logical :: Gra_log,i_print_W,spstmL
    integer,allocatable ::  Q_neigh(:,:)
    integer ::  filen_kt_acc(2)
    
    ! initialize the variables
!    call rw_MC(n_Tsteps,n_sizerelax,n_thousand,restart_MC_steps,Total_MC_Steps,T_relax,T_auto,cone,)
    call rw_MC(io_MC)
    cone=io_MC%cone
    size_table=io_MC%n_Tsteps

    write(6,'(/,a,I6,a,/)') "you are calculating",size_table," temperatures"
    
    allocate(E_av(size_table))
    !     Errors for energy and magnetization
    allocate(E_err_av(size_table),M_err_av(3,size_table))
    !     Energysum and Magnetizationsum and their squaresum
    !     specific heat and suszeptibility
    allocate(C_av(size_table),chi_Q(4,size_table))
    ! everything for the topological charge
    allocate(Q_sq_sum_av(size_table),Qp_sq_sum_av(size_table),Qm_sq_sum_av(size_table))
    !     magnetisation
    allocate(M_sum_av(3,size_table),M_sq_sum_av(3,size_table),chi_l(3,size_table),chi_M(3,size_table),M_av(3,size_table))
    allocate(vortex_av(3,size_table),qeulerp_av(size_table),qeulerm_av(size_table))
    allocate(kt_all(size_table),E_sum_av(size_table),E_sq_sum_av(size_table))
    
    ! updating data during the simulation
    !shape_spin=shape(my_lattices%ordpar%l_modes)
    qeulerp=0.0d0
    qeulerm=0.0d0
    vortex=0.0d0
    magnetization=0.0d0
    E_total=0.0d0
    E_decompose=0.0d0
    kt_all=0.0d0
    Gra_log=io_simu%io_Xstruct
    i_print_W=io_simu%io_warning
    kTini=ext_param%ktini%value
    kTfin=ext_param%ktfin%value
    N_cell=my_lattice%Ncell
    spstmL=io_simu%io_spstmL
    filen_kt_acc=5
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!! Prepare Q_neigh for topological neighbor calculation
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    Call neighbor_Q(my_lattice,Q_neigh)
    E_total=energy_all(Hams,my_lattice)
    call CalculateAverages(my_lattice,Q_neigh,qeulerp,qeulerm,vortex,Magnetization)
    
    ! Measured data
    E_av=0.0d0
    C_av=0.0d0
    chi_Q=0.0d0
    E_err_av=0.0d0
    E_sum_av=0.0d0
    M_err_av=0.0d0
    M_sum_av=0.0d0
    vortex_av=0.0d0
    qeulerp_av=0.0d0
    qeulerm_av=0.0d0
    M_sq_sum_av=0.0d0
    E_sq_sum_av=0.0d0
    Q_sq_sum_av=0.0d0
    Qp_sq_sum_av=0.0d0
    Qm_sq_sum_av=0.0d0
    chi_l=0.0d0
    chi_M=0.0d0
    M_av=0.0d0
    
    ! statistics on the MC
    acc=0.0d0
    rate=0.0d0
    tries=0.0d0
    nb=0.0d0
    
    ! initialize the temperatures
    call ini_temp(kt_all,kTfin,kTini,size_table,i_print_W)
    filen_kt_acc(1)=max(int(log10(maxval(kt_all))),1)
    
    
    write(6,'(/a)') 'Initialization done'
    write(6,'((a,3x,f14.5,a,3x,f9.5))') ' Total energy ', E_total, 'eV    Topological charge ', (qeulerp+qeulerm)/pi*0.25d0
    write(6,'(a/)') '---------------------'
    
    Do n_kT=1,io_MC%n_Tsteps
      kt=kt_all(n_kT)
      call Relaxation(my_lattice,io_MC,N_cell,E_total,E_decompose,Magnetization,qeulerp,qeulerm,kt,acc,rate,nb,cone,hams,Q_neigh)
    !Monte Carlo steps, calculate the values
        do n_MC=1+io_MC%restart_MC_steps,io_MC%Total_MC_Steps+io_MC%restart_MC_steps
    !      Monte Carlo steps for independency
           Do i_relax=1,io_MC%T_auto*N_cell
              Call MCstep(my_lattice,io_MC,N_cell,E_total,E_decompose,Magnetization,kt,acc,rate,nb,cone,Hams)
           End do
           Call MCstep(my_lattice,io_MC,N_cell,E_total,E_decompose,Magnetization,kt,acc,rate,nb,cone,Hams)
    
    ! CalculateAverages makes the averages from the sums
           Call CalculateAverages(my_lattice,Q_neigh,qeulerp_av(n_kT),qeulerm_av(n_kT),Q_sq_sum_av(n_kT),Qp_sq_sum_av(n_kT),Qm_sq_sum_av(n_kT),vortex_av(:,n_kT),vortex &
                    &  ,E_sum_av(n_kT),E_sq_sum_av(n_kT),M_sum_av(:,n_kT),M_sq_sum_av(:,n_kT),E_total,Magnetization)
    
    !**************************
    !!!!!!!!!!!!!!!!!!!!!!!!!!!
       end do ! over n_MC
    !***************************
    !!!!!!!!!!!!!***************************!!!!!!
    !!!!!!!!!!!!!***************************!!!!!!
    
       if (io_MC%n_Tsteps.ne.0) call Calculate_thermo(io_MC%Cor_log,io_MC%total_MC_steps, &
        &    dble(N_cell),kT_all(n_kt),E_sq_sum_av(n_kt),E_sum_av(n_kt),M_sq_sum_av(:,n_kt), &
        &    C_av(n_kt),chi_M(:,n_kt),E_av(n_kt),E_err_av(n_kt),M_err_av(:,n_kt),qeulerp_av(n_kt),qeulerm_av(n_kt),vortex_av(:,n_kt), &
        &    Q_sq_sum_av(n_kt),Qp_sq_sum_av(n_kt),Qm_sq_sum_av(n_kt), &
        &    chi_Q(:,n_kt),chi_l(:,n_kt), &
        &    M_sum_av(:,n_kt),M_av(:,n_kt))
    
    
       write(6,'(5(a,f18.9,2x))') 'M= ',norm(M_av(:,n_kt)), &
         & 'E= ',E_av(n_kT),'Q+= ',qeulerp_av(n_kT),'Q-= ',qeulerm_av(n_kT),'Q= ',qeulerp_av(n_kT)+qeulerm_av(n_kT)
    
    
       if (Gra_log) then
          Call write_config('MC',kt,my_lattice,filen_kt_acc)
       endif
    
      Write(6,'(I6,a,I6,a,f8.4,a,/)')  n_kT, ' nd step out of ', io_MC%n_Tsteps,' steps. T=', kT/k_B,' Kelvin'
    
    end do !over n_kT
    !--------------------------------------------------------------
    !!!!!!!!!!!!!***************************!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!***************************!!!!!!!!!!!!!!!!!!!!!!!

    ! The program of Tobias is used only at last iteration
    if (spstmL) call spstm

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!! write EM.dat output
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    OPEN(7,FILE='EM.dat',action='write',status='unknown', &
      & position='append',form='formatted')
      Write(7,'(27(a,15x))') '# 1:T','2:E_av','3:E_err','4:C','5:M','6:Mx','7:My','8:Mz', &
      & '9:M_err_x','10:M_err_y','11:M_err_z','12:chi_x','13:chi_y','14:chi_z','15:vx', &
      & '16:vy','17:vz','18:qeuler','19:Chi_q','20:Q+','21:Chi_qp','22:Q-','23:Chi_qm', &
      & '24:l_x','25:l_y','26:l_z','27:Chi_QpQm'
    ! write the data into a file
    do i=1,io_MC%n_Tsteps
       Write(7,'(27(E20.10E3,2x),E20.10E3)') kT_all(i)/k_B ,E_av(i), E_err_av(i), C_av(i), norm(M_sum_av(:,i))/N_cell/(io_MC%Total_MC_Steps+io_MC%restart_MC_steps), &
         &             M_sum_av(:,i)/N_cell/(io_MC%Total_MC_Steps+io_MC%restart_MC_steps), M_err_av(:,i), chi_M(:,i), vortex_av(:,i), qeulerm_av(i)+qeulerp_av(i), chi_Q(1,i), &
         &             qeulerp_av(i), chi_Q(2,i), qeulerm_av(i), chi_Q(3,i), chi_l(:,i), chi_Q(4,i)
    enddo
    close(7) 

end subroutine montecarlo
end module
