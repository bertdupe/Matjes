module m_solver_commun
use m_solver
use m_info_dynamics
use m_propagator
use m_eval_BTeff
use m_update_time
use m_solver_order

abstract interface

    function propagator_old(B,damping,mode,size_B)result(propagator)
      integer, intent(in) :: size_B
      real(kind=8), intent(in) :: B(:),damping,mode(:)
      real(kind=8) :: propagator(size_B)
    end function 

    subroutine propagator(B,damping,M,Mout)
        real(8),intent(in)                          ::  damping
        real(8),intent(in),target,contiguous        ::  M(:,:),B(:,:)
        real(8),intent(inout),target,contiguous     ::  Mout(:,:)
    end subroutine

    function integrator_old(Mode_t,D_Mode,BT,dt,size_mode) result(integrator)
      integer, intent(in) :: size_mode
      real(kind=8), intent(in) :: Mode_t(:),D_Mode(:),BT(:),dt
      real(kind=8) :: integrator(size_mode)
    end function 

    function integrator(m,Dmag_int,dt)result(Mout)
        real(8),intent(in),target,contiguous    ::  M(:,:),Dmag_int(:,:)
        real(8),intent(in)                      ::  dt
        real(8),target                          ::  Mout(size(M,1),size(M,2))
    end function


    subroutine temperature_field(kt,damping,mode,BT,DT,size_BT)
      integer, intent(in) :: size_BT
      real(kind=8), intent(in) :: kt,mode(:),damping
      real(kind=8), intent(inout) :: BT(:)
      real(kind=8), intent(inout) :: DT(:)
    end subroutine temperature_field

end interface

procedure (propagator_old), pointer, protected :: get_propagator_field_old=>null()
procedure (propagator), pointer, protected :: get_propagator_field=>null()
procedure (integrator), pointer, protected :: get_integrator_field=>null()
procedure (integrator_old), pointer, protected :: get_integrator_field_old=>null()
procedure (temperature_field), pointer, protected :: get_temperature_field=>null()

private
public :: select_propagator,get_propagator_field,get_propagator_field_old,get_integrator_field,get_integrator_field_old,get_temperature_field

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

call choice_solver(integtype)

! different integration types
!-----------------------------------------------
! Euler integration scheme
!-----------------------------------------------

select case (integtype)
  case (1)

    get_propagator_field => LLG
    get_propagator_field_old => LLG_old
    get_integrator_field => euler
    get_integrator_field_old => euler_old

    N_loop=1

    if (kt.gt.1.0d-7) get_temperature_field => langevin_bath

    call get_butcher_explicit(N_loop)

!
!-----------------------------------------------
! Heun integration scheme
!-----------------------------------------------
  case (2)

    get_propagator_field => LLG
    get_propagator_field_old => LLG_old
    get_integrator_field => euler
    get_integrator_field_old => euler_old

    N_loop=2

    if (kt.gt.1.0d-7) get_temperature_field => langevin_bath

    call get_butcher_explicit(N_loop)

!
!-----------------------------------------------
! SIB(implicit solver)
!-----------------------------------------------
!  case (3)
!    stop 'Does not work so far'
!    get_propagator_field => LLG_B
!    get_integrator_field => implicite
!
!    N_loop=2
!
!    if (kt.gt.1.0d-7) get_temperature_field => wiener_bath
!
!    call get_parameter(io,'input','N_order',N_loop)
!    call get_butcher_implicit(N_loop)
!
!-----------------------------------------------
! Runge Kutta order N integration scheme
!-----------------------------------------------
  case (4)
    get_propagator_field => LLG
    get_propagator_field_old => LLG_old
    get_integrator_field => euler
    get_integrator_field_old => euler_old

    N_loop=4
    call get_parameter(io,'input','N_order',N_loop)

    call get_butcher_explicit(N_loop)

!
!-----------------------------------------------
! SIB without temperature and with error control
!-----------------------------------------------
!  case (6)
!    get_propagator_field => LLG_B
!    get_integrator_field => implicite
!
!    N_loop=2
!
!    if (kt.gt.1.0d-7) get_temperature_field => langevin_bath
!
!    call get_butcher_explicit(N_loop)

!
!-----------------------------------------------
! other cases
!-----------------------------------------------
  case default

   stop 'integration type not implemented'

end select

call close_file('input',io)

end subroutine select_propagator

end module m_solver_commun
