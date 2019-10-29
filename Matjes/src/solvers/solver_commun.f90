module m_solver_commun
use m_solver
use m_info_dynamics
use m_propagator
use m_eval_BTeff
use m_update_time

abstract interface

    function propagator(B,damping,mode,size_B)
      integer, intent(in) :: size_B
      real(kind=8), intent(in) :: B(:),damping,mode(:)
      real(kind=8) :: propagator(size_B)
    end function propagator

    function integrator(Mode_t,D_Mode,BT,dt,size_mode)
      integer, intent(in) :: size_mode
      real(kind=8), intent(in) :: Mode_t(:),D_Mode(:),BT(:),dt
      real(kind=8) :: integrator(size_mode)
    end function integrator

    subroutine temperature_field(kt,damping,mode,BT,DT,size_BT)
      integer, intent(in) :: size_BT
      real(kind=8), intent(in) :: kt,mode(:),damping
      real(kind=8), intent(inout) :: BT(:)
      real(kind=8), intent(inout) :: DT(:)
    end subroutine temperature_field

    real(kind=8) function multiply_dt(dt)
      real(kind=8), intent(in) :: dt
    end function multiply_dt

end interface

procedure (propagator), pointer, protected :: get_propagator_field=>null()
procedure (integrator), pointer, protected :: get_integrator_field=>null()
procedure (temperature_field), pointer, protected :: get_temperature_field=>null()
procedure (multiply_dt), pointer, protected :: multiply=>null()

private
public :: select_propagator,get_propagator_field,get_integrator_field,get_temperature_field,multiply

contains

subroutine select_propagator(kt,N_loop)
use m_io_files_utils
use m_io_utils
implicit none
real(kind=8), intent(in) :: kt
integer, intent(out) :: N_loop
!real(kind=8), intent(inout) :: w(:)
! internal variables
integer :: integtype,io

N_loop=0
io=open_file_read('input')
call get_parameter(io,'input','integration',integtype)
call close_file('input',io)

call choice_solver(integtype)

! different integration types
!-----------------------------------------------
! Euler integration scheme
!-----------------------------------------------
select case (integtype)
  case (1)

    get_propagator_field => LLG
    get_integrator_field => euler
    multiply=>multiply_1

    N_loop=1

    if (kt.gt.1.0d-4) get_temperature_field => langevin_bath

!       if (kt.gt.1.0d-10) call calculate_BTeff(stmtemp,kt,BT_point(iomp)%w)

!
!-----------------------------------------------
! Heun integration scheme
!-----------------------------------------------
  case (2)

    get_propagator_field => LLG
    get_integrator_field => euler
    multiply=>multiply_2

    N_loop=2

    if (kt.gt.1.0d-4) get_temperature_field => langevin_bath

!
!-----------------------------------------------
! SIB(implicit solver)
!-----------------------------------------------
  case (3)
    get_propagator_field => LLG_B
    get_integrator_field => implicite
    multiply=>multiply_1

    N_loop=2

    if (kt.gt.1.0d-4) get_temperature_field => wiener_bath
!
!-----------------------------------------------
! SIA and IMP integration scheme
!-----------------------------------------------
!       case (2)
!       call calculate_Beff(iomp,Beff,spin1,h_int,Hamiltonian)
!
!        spin2(iomp)%w=(integrate(timestep_int,spin1(:,i_x,i_y,i_z,i_m),Beff,kt,damping &
!     & ,stmtemp,i_torque,stmtorque,torque_FL,torque_AFL,adia,nonadia,storque,maxh,Ipol,i_x,i_y,i_z,i_m,spin)+ &
!     & spinini(:,i_x,i_y,i_z,i_m))/2.0d0

!
!-----------------------------------------------
! SIA and IMP integration scheme
!-----------------------------------------------
!       case (4)
!       call calculate_Beff(i_x,i_y,i_z,i_m,Beff,spin,shape_spin,mag_lattice,h_int,Hamiltonian)

!        spinafter(:,i_x,i_y,i_z,i_m)=(integrate(timestep_int,spin(4:6,i_x,i_y,i_z,i_m),Beff,kt,damping &
!     & ,stmtemp,i_torque,stmtorque,torque_FL,torque_AFL,adia,nonadia,storque,maxh,check,Ipol,i_x,i_y,i_z,i_m,spin)+ &
!     &  spinini(:,i_x,i_y,i_z,i_m))/2.0d0

!
!-----------------------------------------------
! SIB without temperature and with error control
!-----------------------------------------------
  case (6)
    get_propagator_field => LLG_B
    get_integrator_field => implicite
    multiply=>multiply_1

    N_loop=2

    if (kt.gt.1.0d-4) get_temperature_field => langevin_bath

!
!-----------------------------------------------
! other cases
!-----------------------------------------------
  case default

   stop 'not implemented'

end select

end subroutine select_propagator

end module m_solver_commun
