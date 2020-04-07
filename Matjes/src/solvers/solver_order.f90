module m_solver_order

real(kind=8), allocatable :: Butcher_table(:,:)

private
public :: get_butcher,get_D_mode

contains

!
! subroutine that getsD_mode as a function of the Butcher's table
!
subroutine get_D_mode(D_mode,j,N_loop,D_out)
use m_convert
implicit none
integer, intent(in) :: j,N_loop
real(kind=8), intent(in) :: D_mode(:,:)
real(kind=8), intent(inout) :: D_out(:)
! internal
integer :: i

D_out=0.0d0

do i=1,N_loop
  D_out=D_out+Butcher_table(i,j+1)*D_mode(:,i)
enddo

end subroutine get_D_mode

!
! subroutine that gets the N spliting the the segments of size dt
!
! Butcher 1987, Lambert 1991
subroutine get_butcher(N)
use m_convert
implicit none
integer, intent(in) :: N
! internal
character(len=50) :: form
integer :: i,j

form=convert('(',N,'(2x,f6.4))')

allocate(Butcher_table(N,N+1))
Butcher_table=0.0d0

select case (N)
  case (1)
    write(6,'(a)') 'Euler solver'
    Butcher_table(1,2)=1.0d0

  case (2)
    write(6,'(a)') 'Heun solver'
    Butcher_table(1,2)=1.0d0
    Butcher_table(:,3)=0.5d0

  case (3)
    write(6,'(a)') 'standard Runge Kutta 3'

    Butcher_table(1,2)=0.5d0
    Butcher_table(:,3)=(/1.0d0,0.5d0,0.0d0/)
    Butcher_table(:,4)=(/1.0d0,2.0d0,1.0d0/)/6.0d0

  case (4)
    write(6,'(a)') 'standard Runge Kutta 4'

    Butcher_table(1,2)=0.5d0
    Butcher_table(2,3)=0.5d0
    Butcher_table(3,4)=1.0d0
    Butcher_table(:,5)=(/1.0d0,2.0d0,2.0d0,1.0d0/)/6.0d0

  case (6)
    write(6,'(a)') 'Fehlberg method (order 6)'

    Butcher_table(1,2)=0.25d0
    Butcher_table(1:2,3)=(/3.0d0,9.0d0/)/32.0d0
    Butcher_table(1:3,4)=(/1932.0d0,-7200.0d0,7296.0d0/)/2197.0d0
    Butcher_table(1:4,5)=(/439.0d0/216.0d0,-8.0d0,3680.0d0/513.0d0,-845.0d0/4104.0d0/)
    Butcher_table(1:5,6)=(/-8.0d0/27.0d0,2.0d0,-3544.0d0/2565.0d0,1859.0d0/4104.0d0,-11.0d0/40.0d0/)
    Butcher_table(:,7)=(/25.0d0/216.0d0,0.0d0,1408.0d0/2565.0d0,2197.0d0/4104.0d0,-1.0d0/5.0d0,0.0d0/)

  case (7)
    write(6,'(a)') 'Dormand-Prince method (order 7)'

    Butcher_table(1,2)=1.0d0/5.0d0
    Butcher_table(1:2,3)=(/3.0d0 , 9.0d0/) /40.0d0
    Butcher_table(1:3,4)=(/44.0d0/45.0d0 ,-56.0d0/15.0d0 , 32.0d0/9.0d0/)
    Butcher_table(1:4,5)=(/19372.0d0/6561.0d0 , -25360.0d0/2187.0d0 , 64448.0d0/6561.0d0 , -212.0d0/729.0d0/)
    Butcher_table(1:5,6)=(/9017.0d0/3168.0d0 , -355.0d0/33.0d0 , 46732.0d0/5247.0d0 , 49.0d0/176.0d0 , -5103.0d0/18656.0d0/)
    Butcher_table(1:6,7)=(/35.0d0/384.0d0 , 0.0d0 , 500.0d0/1113.0d0 , 125.0d0/192.0d0 , -2187.0d0/6784.0d0 , 11.0d0/84.0d0/)
    Butcher_table(:,8)=(/5179.0d0/57600.0d0 , 0.0d0 , 7571.0d0/16695.0d0 , 393.0d0/640.0d0 , -92097.0d0/339200.0d0 , 187.0d0/2100.0d0 , 1.0d0/40.0d0/)

  case default
    stop 'get_D_mode - case not coded'

end select


write(6,'(a)') ''
write(6,'(a)') ' Butchers table'
do i=1,N+1
  write(6,form) (Butcher_table(j,i),j=1,N)
enddo
write(6,'(a)') ''

end subroutine get_butcher



end module m_solver_order
