module m_minimize
use m_H_public
implicit none

real(kind=8) :: masse=1.0d0
real(kind=8) :: dt=0.1d0
integer :: N_minimization=1000 !OBSOLETE
real(kind=8) :: conv_torque=1.0d-6

interface minimize
  module procedure minimize_lattice
end interface

interface minimize_infdamp
  module procedure minimize_infdamp_lattice
  !module procedure minimize_infdamp_2Dlattice,minimize_infdamp_lattice
end interface

private
public :: minimize,minimize_infdamp

contains

!
!
! Minimization routine that it actually used.
! The interface is used to put the data into the good format
!
!
subroutine minimize_lattice(lat,io_simu,Hams)
    use m_derived_types, only : io_parameter,lattice
    use m_basic_types, only : vec_point
    use m_constants, only : pi
    use m_write_spin
    use m_createspinfile
    use m_solver, only : minimization
    use m_vector, only : norm_cross,norm, calculate_damping
    use m_dyna_utils, only : copy_lattice
    use m_eval_Beff
    use m_lattice, only : my_order_parameters
    use m_operator_pointer_utils
    use m_Beff_H, only: get_B
    
    implicit none
    type(io_parameter), intent(in) :: io_simu
    type(lattice), intent(inout) :: lat
    class(t_H), intent(in) :: Hams(:)
    ! dummy variable
    real(kind=8),allocatable, dimension(:,:) :: velocity,predicator,force
    real(kind=8),allocatable, dimension(:) :: V_eff,F_temp,Beff
    ! internal
    real(kind=8) :: dumy,force_norm,Energy,vmax,vtest,Eint,test_torque,max_torque
    ! the computation time
    integer :: i_min,gra_freq
    logical :: gra_log,i_magnetic
    integer :: iomp,dim_mode,N_cell

    ERROR STOP "MINIMIZE HAS NOT BEEN UPDATED IN A LONG TIME AND IS PROBABLY OBSOLETE (USE MINIMIZE_INFDAMP?)"
    
    call init_variables()
    gra_freq=io_simu%io_frequency
    gra_log=io_simu%io_Xstruct
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    force_norm=0.0d0
    test_torque=0.0d0
    max_torque=10.0d0
    vmax=0.0d0
    vtest=0.0d0
    N_cell=lat%Ncell
    dim_mode=lat%M%dim_mode
    
    allocate(Beff(dim_mode*N_cell),source=0.0d0)
    allocate(velocity(dim_mode,N_cell),predicator(dim_mode,N_cell),force(dim_mode,N_cell),source=0.0d0)
    allocate(V_eff(dim_mode),F_temp(dim_mode),source=0.0d0)
    !ANGLEICHEN AN MINIMIZE_INFDAMP
    
    Call get_B(Hams,lat,Beff)
!    do iomp=1,N_cell
!    do iomp=1,1
!       call calculate_Beff(F_eff,iomp,lat%ordpar%all_l_modes)
!       force(:,iomp)=calculate_damping(lat%ordpar%all_l_modes(iomp)%w,F_eff)
!       call minimization(lat%ordpar%all_l_modes(iomp)%w,force(:,iomp),predicator(:,iomp),dt**2,masse*2.0d0)
!       test_torque=norm_cross(predicator(:,iomp),force(:,iomp),1,3)
!    enddo
    write(*,*) Beff(1:3)
    STOP "DERP"
#if 0
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
    
      Energy=energy_all(Hams,lat)
    
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
    
    ERROR STOP "FINISH minimize"
    
    nullify(my_lattice,all_mode)
    
#endif
end subroutine


subroutine minimize_infdamp_lattice(lat,io_simu,Hams)
    use m_derived_types, only : io_parameter,lattice
    use m_basic_types, only : vec_point
    use m_constants, only : pi
    use m_write_spin
    use m_createspinfile
    use m_vector, only : cross,norm
    use m_eval_Beff
    use m_lattice, only : my_order_parameters
    use m_operator_pointer_utils
    use m_Beff_H, only: get_B
    type(io_parameter), intent(in)  :: io_simu
    type(lattice),intent(inout)     :: lat
    class(t_H), intent(in)          :: Hams(:)
    ! internal
    real(8)                     :: max_torque,test_torque,Edy
    integer                     :: iomp,iter,gra_freq
    logical                     :: gra_log
    integer                     :: N_cell,N_dim,N_mag
    real(8),allocatable,target  :: F_eff(:)
    real(8),allocatable         :: F_norm(:),torque(:)
    real(8),pointer             :: M3(:,:),F_eff3(:,:)
    
    write(6,'(/,a,/)') 'entering the infinite damping minimization routine'
    
    call init_variables()
    N_cell=lat%Ncell
    N_dim=lat%M%dim_mode
    N_mag=(N_dim/3)*N_cell
    
    allocate(F_norm(N_mag),source=0.0d0)
    allocate(F_eff(N_dim*N_cell),torque(N_dim*N_cell),source=0.0d0)
    F_eff3(1:3,1:N_mag)=>F_eff
    M3(1:3,1:N_mag)=>lat%M%all_modes
    
    gra_log=io_simu%io_Xstruct
    gra_freq=io_simu%io_frequency
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Prepare the calculation of the energy and the effective field
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    Edy=energy_all(Hams,lat)
    write(6,'(/a,2x,E20.12E3/)') 'Initial total energy density (eV/fu)',Edy/real(N_cell,8)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !           Begin minimization
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    iter=0
    max_torque=10.0d0
    do while (max_torque.gt.conv_torque)
        max_torque=0.0d0

        Call get_B(Hams,lat,F_eff)
        F_norm=norm2(F_eff3,1)
        if(any(F_norm.lt.1.0d-8)) stop 'problem in the infinite damping minimization routine' !avoid divide by 0
        torque=cross(lat%M%all_modes,F_eff,1,N_dim*N_cell)
        test_torque=max(maxval(torque),-minval(torque))
        if ( abs(test_torque).gt.max_torque ) max_torque=test_torque
        !align the moments onto normalized field
        do iomp=1,N_mag
            M3(:,iomp)=F_eff3(:,iomp)/F_norm(iomp)
        enddo
        iter=iter+1
        !print max_torque every 100 iterations
        if (mod(iter,100).eq.0) write(*,*) 'Max torque =',max_torque
    
        !write config to files
        if ((gra_log).and.(mod(iter,gra_freq).eq.0)) then
             call CreateSpinFile(iter/gra_freq,lat%ordpar%all_l_modes)
             call WriteSpinAndCorrFile(iter/gra_freq,lat%ordpar%all_l_modes,'SpinSTM_')
             write(6,'(a,3x,I10)') 'wrote Spin configuration and povray file number',iter/gra_freq
        endif
    enddo
    write(*,*) 'Max_torque=',max_torque,' tolerance reached, minimization completed in ',iter,' iterations.'
    Edy=energy_all(Hams,lat)
    write(6,'(/a,2x,E20.12E3/)') 'Final total energy density (eV/fu)',Edy/real(N_cell,8)
    nullify(M3,F_eff3)
    deallocate(F_eff,F_norm,torque)
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
