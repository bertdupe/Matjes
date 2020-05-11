module m_spline

contains

subroutine spline_hermite_set ( ndata, tdata, c )
!

  implicit none
!
  integer, intent(in) :: ndata
!
  real(kind=8), intent(inout) :: c(4,ndata)
  real(kind=8) :: divdif1, divdif3,dt
  integer i
  real(kind=8), intent(in) :: tdata(ndata)
!
  do i = 1, ndata-1
    dt = tdata(i+1) - tdata(i)
    divdif1 = ( c(1,i+1) - c(1,i) ) / dt
    divdif3 = c(2,i) + c(2,i+1) - 2d0* divdif1
    c(3,i) = ( divdif1 - c(2,i) - divdif3 ) / dt
    c(4,i) = divdif3 / (dt*dt)
  end do

  c(3,ndata) = 0d0
  c(4,ndata) = 0d0


end subroutine spline_hermite_set


subroutine spline_hermite_val ( ndata, tdata, c, tval, sval, dsval )
!

!
  implicit none
!
  integer, intent(in) :: ndata
!
  real(kind=8), intent(in) :: c(4,ndata),tdata(ndata),tval
  real(kind=8) :: dt
  integer :: left
  real(kind=8), intent(out) :: sval, dsval

!
!  Find the interval [ TDATA(LEFT), TDATA(RIGHT) ] that contains
!  or is nearest to TVAL.
!
  call rvec_bracket ( ndata, tdata, tval, left)
!
!  Evaluate the cubic polynomial.
!
  dt = tval - tdata(left)

  sval = c(1,left) + dt * ( c(2,left) + dt * ( c(3,left) + dt * c(4,left) ) )
  dsval = c(2,left) + dt*(2d0*c(3,left) + dt*3d0*c(4,left))


end subroutine spline_hermite_val





subroutine rvec_bracket ( n, x, xval, left)
!
implicit none
!
integer, intent(in) :: n
real(kind=8), intent(in) :: x(n),xval
integer, intent(out) :: left
! internal variables
integer i

!
do i = 2, n - 1

   if ( xval < x(i) ) then
      left = i - 1

      return
   end if

end do

left = n - 1

end subroutine rvec_bracket





subroutine hermite_fit(n,nn,x,y,dydx,xx,yy,dyy,c)
implicit none
integer, intent(in) :: n,nn
real(kind=8), intent(in) :: x(n),y(n),dydx(n)
real(kind=8), intent(out) :: xx(nn),yy(nn),dyy(nn),c(4,n)
! internal variables
real(kind=8) :: dx
integer :: i

do i =1,n
   c(1,i) = y(i)
   c(2,i) = dydx(i)
end do

dx = (x(n)-x(1))/(nn-1)

call spline_hermite_set(n,x,c)

do i=1,nn
   xx(i) = x(1) + (i-1)*dx
   call spline_hermite_val(n,x,c,xx(i),yy(i),dyy(i))
end do

end subroutine hermite_fit



end module m_spline
