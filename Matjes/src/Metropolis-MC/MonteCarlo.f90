module m_montecarlo
implicit none
contains
!
! Routine that does the Monte Carlo (and not the parallel tempering)
!
!subroutine MonteCarlo(my_lattices,mag_motif,io_simu,gra_topo,ext_param)
subroutine montecarlo(my_lattice,motif,io_simu,ext_param,Hams)
use m_constants, only : k_b,pi
use m_vector, only : norm
use m_derived_types, only : lattice,t_cell,io_parameter,simulation_parameters,point_shell_Operator
use m_modes_variables, only : point_shell_mode
use m_basic_types, only : vec_point
use m_rw_MC
use m_lattice, only : my_order_parameters
use m_topo_commons
use m_convert
use m_io_files_utils
use m_operator_pointer_utils
use m_eval_Beff
use m_MCstep
use m_Htype_gen
use m_relaxation
!use m_H_type

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
#ifdef CPP_MPI
      use m_mpi_prop, only : MPI_COMM
      use m_gather_reduce
#endif
type(lattice), intent(inout) :: my_lattice
type(t_cell), intent(in) :: motif
type(io_parameter), intent(in) :: io_simu
type(simulation_parameters), intent(in) :: ext_param
class(t_H), intent(in) :: Hams(:)

!!!!!!!!!!!!!!!!!!!!!!!
! internal variables
!!!!!!!!!!!!!!!!!!!!!!!
! slope of the MC
integer :: i_relax,n_kT,n_MC,T_relax,restart_MC_steps,T_auto,N_cell,n_sizerelax,n_thousand,n_Tsteps,size_table,Total_MC_Steps
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
logical :: underrel,overrel,sphere,equi,i_restart,ising,print_relax,Cor_log,Gra_log,i_magnetic,i_print_W,spstmL,i_temperature

type(vec_point),pointer :: all_mode(:)
type(vec_point),allocatable,dimension(:) :: mode_magnetic,mode_temp

! initialize the variables
call rw_MC(n_Tsteps,n_sizerelax,n_thousand,restart_MC_steps,Total_MC_Steps,T_relax,T_auto,cone,i_restart,ising,underrel,overrel,sphere,equi,print_relax,Cor_log)
size_table=n_Tsteps

#ifdef CPP_MPI
if (i_separate) size_table=isize*nRepProc
!      if (i_ghost) size_table=isize/n_ghost
if (irank_working.eq.0) write(6,'(/,a,I6,a,/)') "you are calculating",size_table," temperatures"
#else
size_table=n_Tsteps
write(6,'(/,a,I6,a,/)') "you are calculating",size_table," temperatures"
#endif

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
N_cell=product(shape(my_lattice%ordpar%l_modes))
spstmL=io_simu%io_spstmL

all_mode=>my_lattice%ordpar%all_l_modes
!
! Prepare the calculation of the Energy and the field
! magnetization
do i=1,size(my_order_parameters)
  if ('magnetic'.eq.trim(my_order_parameters(i)%name)) then
   allocate(mode_magnetic(N_cell))
   call dissociate(mode_magnetic,N_cell)
   call associate_pointer(mode_magnetic,all_mode,'magnetic',i_magnetic)
  endif
enddo

! temperature
do i=1,size(my_order_parameters)
  if ('temperature'.eq.trim(my_order_parameters(i)%name)) then
   allocate(mode_temp(N_cell))
   call dissociate(mode_temp,N_cell)
   call associate_pointer(mode_temp,all_mode,'temperature',i_temperature)

   exit
  endif
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! associate pointer for the topological charge, vorticity calculations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call get_size_Q_operator(my_lattice)
call associate_Q_operator(all_mode,my_lattice%boundary,shape(my_lattice%ordpar%l_modes))

!ktini=ext_param%ktini%value
!ktfin=ext_param%ktfin%value
! initializing the variables above

E_total=energy_all(Hams,my_lattice)

call CalculateAverages(mode_magnetic,qeulerp,qeulerm,vortex,Magnetization)

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
#ifdef CPP_MPI
call ini_temp(kt_all,kTfin,kTini,size_table,irank_working,nRepProc,i_print_W)
#else
call ini_temp(kt_all,kTfin,kTini,size_table,i_print_W)
#endif


write(6,'(/a)') 'Initialization done'
write(6,'(2(a,3x,f9.5))') ' Total energy ', E_total, 'eV    Topological charge ', qeulerp+qeulerm
write(6,'(a/)') '---------------------'

Do n_kT=1,n_Tsteps

  kt=kt_all(n_kT)
!  do i=1,N_cell
!     mode_temp(i)%w=kt
!  enddo
!  call set_E_matrix(my_lattice%dim_mode,kt/k_b)

  call Relaxation(my_lattice,N_cell,n_sizerelax,n_thousand,T_relax,E_total,E_decompose,Magnetization,qeulerp,qeulerm,kt,acc,rate,nb,cone,ising,equi,overrel,sphere,underrel,print_relax,hams)

!       Monte Carlo steps, calculate the values

    do n_MC=1+restart_MC_steps,Total_MC_Steps+restart_MC_steps

!         Monte Carlo steps for independency

       Do i_relax=1,T_auto*N_cell

          Call MCstep(my_lattice,N_cell,E_total,E_decompose,Magnetization,kt,acc,rate,nb,cone,ising,equi,overrel,sphere,underrel,Hams)

       End do

       Call MCstep(my_lattice,N_cell,E_total,E_decompose,Magnetization,kt,acc,rate,nb,cone,ising,equi,overrel,sphere,underrel,Hams)

! Calculate the topological charge and the vorticity
       dumy=get_charge()
       qeulerp_av(n_kT)=qeulerp_av(n_kT)+dumy(1)
       qeulerm_av(n_kT)=qeulerm_av(n_kT)+dumy(2)
       vortex_av(:,n_kT)=vortex_av(:,n_kT)+dumy(3:5)

! CalculateAverages makes the averages from the sums
       Call CalculateAverages(qeulerp_av(n_kT),qeulerm_av(n_kT),Q_sq_sum_av(n_kT),Qp_sq_sum_av(n_kT),Qm_sq_sum_av(n_kT),vortex_av(:,n_kT),vortex &
                &  ,E_sum_av(n_kT),E_sq_sum_av(n_kT),M_sum_av(:,n_kT),M_sq_sum_av(:,n_kT),E_total,Magnetization)
!      if (Cor_log) chi_l(:,n_kT)=chi_l(:,n_kT)+Correlation(spin_sum,spin(4:6,:,:,:,:),shape_spin,n_MC,dble(N_cell))


! Calculate the sum of the spin components and angles for average
!      if (CalTheta) Call SphericalCoordinates(spin(4:6,:,:,:,:),shape_spin,angle_sum)

!**************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!
   end do ! over n_MC
!***************************
!!!!!!!!!!!!!***************************!!!!!!
!!!!!!!!!!!!!***************************!!!!!!

   if (n_Tsteps.ne.0) call Calculate_thermo(Cor_log,total_MC_steps, &
    &    dble(N_cell),kT_all(n_kt),E_sq_sum_av(n_kt),E_sum_av(n_kt),M_sq_sum_av(:,n_kt), &
    &    C_av(n_kt),chi_M(:,n_kt),E_av(n_kt),E_err_av(n_kt),M_err_av(:,n_kt),qeulerp_av(n_kt),qeulerm_av(n_kt),vortex_av(:,n_kt), &
    &    Q_sq_sum_av(n_kt),Qp_sq_sum_av(n_kt),Qm_sq_sum_av(n_kt), &
    &    chi_Q(:,n_kt),chi_l(:,n_kt), &
    &    M_sum_av(:,n_kt),M_av(:,n_kt))


   write(6,'(5(a,f18.9,2x))') 'M= ',norm(M_av(:,n_kt)), &
     & 'E= ',E_av(n_kT),'Q+= ',qeulerp_av(n_kT),'Q-= ',qeulerm_av(n_kT),'Q= ',qeulerp_av(n_kT)+qeulerm_av(n_kT)


   if (Gra_log) then
      call CreateSpinFile(kt,mode_magnetic)
      call WriteSpinAndCorrFile(kt,all_mode,'SpinSTM_')

   endif

!ccccccccccccccccccccccccccccccccccccc
! Calculate the topological charge
!cccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccc
! Calculate the oriented topological charge
!cccccccccccccccccccccccccccccccccccc


! ===========================================================
! Write the averaged coordinate of spin in spherical coordinates
! ===========================================================

! ===========================================================
! Computation of the energy density on the lattice
! ===========================================================

!   if (CalEnergy) then
!
!      write(fname,'(a,14a,a)')'Energy_T_',(toto(i:i),i=1, &
!            &  len_trim(toto)),'.dat'
!      write(fname2,'(a,14a,a)')'DensityOfEnergy_T_',(toto(i:i),i=1, &
!            &  len_trim(toto)),'.dat'
!
!      Call EnergyDensity(fname,fname2,spin,shape_spin,tableNN,shape_tableNN,masque,shape_masque,indexNN,shape_index)
!   endif

  Write(6,'(I6,a,I6,a,f8.4,a,/)')  n_kT, ' nd step out of ', n_Tsteps,' steps. T=', kT/k_B,' Kelvin'

end do !over n_kT
!--------------------------------------------------------------
!!!!!!!!!!!!!***************************!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!***************************!!!!!!!!!!!!!!!!!!!!!!!

! The program of Tobias is used only at last iteration
#ifdef CPP_MPI
        if ((irank_working.eq.0).AND.spstmL) call spstm
#else
        if (spstmL) call spstm
#endif

#ifdef CPP_MPI
    if (irank_working.eq.0) then
#endif
      OPEN(7,FILE='EM.dat',action='write',status='unknown', &
        & position='append',form='formatted')
        Write(7,'(27(a,15x))') '# 1:T','2:E_av','3:E_err','4:C','5:M','6:Mx','7:My','8:Mz', &
        & '9:M_err_x','10:M_err_y','11:M_err_z','12:chi_x','13:chi_y','14:chi_z','15:vx', &
        & '16:vy','17:vz','18:qeuler','19:Chi_q','20:Q+','21:Chi_qp','22:Q-','23:Chi_qm', &
        & '24:l_x','25:l_y','26:l_z','27:Chi_QpQm'
#ifdef CPP_MPI
    endif
#endif

#ifdef CPP_MPI

if (i_separate) call end_gather(kT_all,E_av,E_err_av,C_av,M_sum_av,M_err_av,chi_M,vortex_av,chi_Q,qeulerp_av,qeulerm_av,chi_l, &
                          & size_table,irank_working,n_Tsteps,MPI_COMM)

   ! case of the ghost
!        if (i_ghost) then
!            n_Tsteps=isize/n_ghost
!            if (irank_box.eq.0) call end_gather(shape_masque,vortex_mpi,qeulerp_mpi,qeulerm_mpi,Im_mpi,Qim_mpi,kt_mpi,mpi_l, &
!            &     E_av(:),C(:),M_err(:),E_err(:),chi(:),M(:,:),n_Tsteps,MPI_COMM_MASTER)
!        endif

!        if ((i_separate.or.i_paratemp).and.(.not.i_ghost)) then
!            n_Tsteps=isize

!            call end_gather(shape_masque,vortex_mpi,qeulerp_mpi,qeulerm_mpi,Im_mpi,Qim_mpi,kt_mpi,mpi_l, &
!            &    E_av(:),C(:),M_err(:),E_err(:),chi(:),M(:,:),isize,MPI_COMM)
!        endif

!        if ((irank_working.eq.0).and.(.not.i_paratemp)) then
!            do i=1,n_Tsteps
!                Write(7,'(23(E20.10E3,2x),E20.10E3)') kt_mpi(i)/k_B,E_av(i),E_err(i),C(i), &
!                & norm(M(:,i)), M(:,i),M_err(:,i), chi(i), vortex_mpi(:,i) &
!                & ,qeulerp_mpi(i)+qeulerm_mpi(i),qeulerp_mpi(i),qeulerm_mpi(i), &
!                & Im_mpi(i,1),Qim_mpi(i,1),mpi_l(:,i)
!            enddo
!        endif
#endif

#ifdef CPP_MPI
if (irank_working.eq.0) then

   do i=1,size_table
      Write(7,'(27(E20.10E3,2x),E20.10E3)') kT_all(i)/k_B ,E_av(i), E_err_av(i), C_av(i), norm(M_sum_av(:,i))/N_cell/(Total_MC_Steps+restart_MC_steps), &
     &             M_sum_av(:,i), M_err_av(:,i), chi_M(:,i), vortex_av(:,i), qeulerm_av(i)+qeulerp_av(i), chi_Q(1,i), &
     &             qeulerp_av(i), chi_Q(2,i), qeulerm_av(i), chi_Q(3,i), chi_l(:,i), chi_Q(4,i)
   enddo

   close(7) !Close EM.dat

endif
#else
! write the data into a file
do i=1,n_Tsteps
   Write(7,'(27(E20.10E3,2x),E20.10E3)') kT_all(i)/k_B ,E_av(i), E_err_av(i), C_av(i), norm(M_sum_av(:,i))/N_cell/(Total_MC_Steps+restart_MC_steps), &
     &             M_sum_av(:,i)/N_cell/(Total_MC_Steps+restart_MC_steps), M_err_av(:,i), chi_M(:,i), vortex_av(:,i), qeulerm_av(i)+qeulerp_av(i), chi_Q(1,i), &
     &             qeulerp_av(i), chi_Q(2,i), qeulerm_av(i), chi_Q(3,i), chi_l(:,i), chi_Q(4,i)
enddo

close(7) !Close EM.dat
#endif

nullify(all_mode)
end subroutine montecarlo
end module
