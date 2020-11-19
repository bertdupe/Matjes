module m_solver_order

real(kind=8), allocatable :: Butcher_table(:,:)

private
public :: get_butcher_explicit,get_D_mode,get_dt_mode,get_butcher_implicit
public :: get_Dmag_int

contains

!
! subroutine that checks the new version of dt depending on the Butcher's table
! if dt=0 do not carry out integration
real(kind=8) function get_dt_mode(dt,i_loop)
implicit none
integer :: i_loop
real(kind=8) :: dt
! internal

get_dt_mode=Butcher_table(1,i_loop+1)*dt

end function get_dt_mode

!
! subroutine that getsD_mode as a function of the Butcher's table
!
subroutine get_D_mode(D_mode,j,N_loop,D_out)
implicit none
integer, intent(in) :: j,N_loop
real(kind=8), intent(in) :: D_mode(:,:)
real(kind=8), intent(inout) :: D_out(:)
! internal
integer :: i

D_out=0.0d0

do i=2,N_loop+1
  D_out=D_out+Butcher_table(i,j+1)*D_mode(:,i-1)
enddo

end subroutine get_D_mode

!
! subroutine that gets D_mode as a function of the Butcher's table
!
subroutine get_Dmag_int(Dmag,i_loop,N_loop,Dmag_int)
    real(8),intent(in)      ::  Dmag(:,:,:)
    real(8),intent(out)     ::  Dmag_int(:,:)
    integer,intent(in)      ::  i_loop,N_loop

    integer                 ::  i

    Dmag_int=0.d0
    do i=2,N_loop+1
        Dmag_int=Dmag_int+Butcher_table(i,i_loop+1)*Dmag(:,:,i-1)
    enddo
end subroutine


!
! subroutine that gets the N spliting the the segments of size dt
!
! Butcher 1987, Lambert 1991
subroutine get_butcher_explicit(N)
implicit none
integer, intent(in) :: N
! internal
integer :: i,j

allocate(Butcher_table(N+1,N+1))
Butcher_table=0.0d0

select case (N)
  case (1)
    write(6,'(a)') 'Euler solver'
    Butcher_table(:,2)=1.0d0

  case (2)
    write(6,'(a)') 'Heun solver'
    Butcher_table(1:2,2)=1.0d0
    Butcher_table(:,3)=(/1.0d0,0.5d0,0.5d0/)

  case (3)
    write(6,'(a)') 'standard Runge Kutta 3'

    Butcher_table(1:2,2)=0.5d0
    Butcher_table(:,3)=(/1.0d0,-1.0d0,2.0d0,0.0d0/)
    Butcher_table(:,4)=(/6.0d0,1.0d0,2.0d0,1.0d0/)/6.0d0

  case (4)
    write(6,'(a)') 'standard Runge Kutta 4'

    Butcher_table(1:2,2)=0.5d0
    Butcher_table(:,3)=(/0.5d0,0.0d0,0.5d0,0.0d0,0.0d0/)
    Butcher_table(:,4)=(/1.0d0,0.0d0,0.0d0,1.0d0,0.0d0/)
    Butcher_table(:,5)=(/6.0d0,1.0d0,2.0d0,2.0d0,1.0d0/)/6.0d0

  case (6)
    write(6,'(a)') 'Fehlberg method (order 6)'

    Butcher_table(1:2,2)=0.25d0
    Butcher_table(1:3,3)=(/12.0d0,3.0d0,9.0d0/)/32.0d0
    Butcher_table(1:4,4)=(/2028.0d0,1932.0d0,-7200.0d0,7296.0d0/)/2197.0d0
    Butcher_table(1:5,5)=(/1.0d0,439.0d0/216.0d0,-8.0d0,3680.0d0/513.0d0,-845.0d0/4104.0d0/)
    Butcher_table(1:6,6)=(/0.5d0,-8.0d0/27.0d0,2.0d0,-3544.0d0/2565.0d0,1859.0d0/4104.0d0,-11.0d0/40.0d0/)
    Butcher_table(:,7)=(/1.0d0,25.0d0/216.0d0,0.0d0,1408.0d0/2565.0d0,2197.0d0/4104.0d0,-1.0d0/5.0d0,0.0d0/)

  case (7)
    write(6,'(a)') 'Dormand-Prince method (order 7)'

    Butcher_table(1:2,2)=1.0d0/5.0d0
    Butcher_table(1:3,3)=(/12.0d0, 3.0d0 , 9.0d0/) /40.0d0
    Butcher_table(1:4,4)=(/4.0d0/5.0d0, 44.0d0/45.0d0 ,-56.0d0/15.0d0 , 32.0d0/9.0d0/)
    Butcher_table(1:5,5)=(/8.0d0/9.0d0, 19372.0d0/6561.0d0 , -25360.0d0/2187.0d0 , 64448.0d0/6561.0d0 , -212.0d0/729.0d0/)
    Butcher_table(1:6,6)=(/1.0d0, 9017.0d0/3168.0d0 , -355.0d0/33.0d0 , 46732.0d0/5247.0d0 , 49.0d0/176.0d0 , -5103.0d0/18656.0d0/)
    Butcher_table(1:7,7)=(/1.0d0, 35.0d0/384.0d0 , 0.0d0 , 500.0d0/1113.0d0 , 125.0d0/192.0d0 , -2187.0d0/6784.0d0 , 11.0d0/84.0d0/)
    Butcher_table(:,8)=(/1.0d0, 5179.0d0/57600.0d0 , 0.0d0 , 7571.0d0/16695.0d0 , 393.0d0/640.0d0 , -92097.0d0/339200.0d0 , 187.0d0/2100.0d0 , 1.0d0/40.0d0/)

  case default
    stop 'get_D_mode - case not coded'

end select

call print_Butcher(N+1)

end subroutine get_butcher_explicit


!
!
!
!
subroutine get_butcher_implicit(N)
implicit none
integer, intent(in) :: N
! internal
integer :: i,j

allocate(Butcher_table(N+1,N+1))
Butcher_table=0.0d0

select case (N)
  case (1)
    write(6,'(a)') 'Euler solver'
    Butcher_table=1.0d0

  case default
    stop 'get_D_mode - case not coded'

end select

call print_Butcher(N+1)

end subroutine get_butcher_implicit


!
!
!
!
subroutine print_Butcher(N)
use m_convert
implicit none
integer, intent(in) :: N
! internal
integer :: i,j
character(len=50) :: form

form=convert('(',N+1,'(2x,f10.5))')

write(6,'(a)') ''
write(6,'(a)') ' Butchers table'
do i=1,N
  write(6,form) (Butcher_table(j,i),j=1,N)
enddo
write(6,'(a)') ''

end subroutine print_Butcher

end module m_solver_order
