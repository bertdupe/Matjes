module m_matrix_power
  implicit none

! Overloading the ** operator does not work because the compiler cannot
! differentiate between matrix exponentiation and the elementwise raising
! of an array to a power therefore we define a new operator
  interface operator (.matpow.)
    module procedure matrix_exp
  end interface

contains

function matrix_exp(m, n) result (res)
real(kind=8), intent(in)  :: m(:,:)
integer, intent(in)  :: n
real(kind=8) :: res(size(m,1),size(m,2))
! internal
integer :: i

if(n == 0) then
  res = 0
  do i = 1, size(m,1)
    res(i,i) = 1.0d0
  end do
  return
end if

res = m
do i = 2, n
  res = matmul(res, m)
end do

end function matrix_exp

end module m_matrix_power
