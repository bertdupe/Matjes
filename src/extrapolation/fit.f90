module m_fit
use m_fit_public
use m_fit_base
use m_derived_types, only : lattice
implicit none

private
public :: execute_fit

contains

subroutine execute_fit(pos,my_lattice)
  type(lattice), intent(in) :: my_lattice
  real(8), intent(in) :: pos(:)

  ! internal
  class(fit), allocatable :: my_fit
  real(8),allocatable     :: X(:,:),Y(:)

  integer :: N,N_variable,i

  N=size(pos)/3
  allocate(X(2,N),source=0.0d0)
  allocate(Y(N),source=0.0d0)

  do i=1,N
     X(1,i)=pos(3*(i-1)+1)
     X(2,i)=pos(3*(i-1)+2)
     Y(i)=my_lattice%M%modes_3(3,i)
  enddo

  call choose_fit(my_fit)
  call my_fit%init()
  call my_fit%execute(X,Y,0.0d0)
  call my_fit%init_base()

end subroutine

end module
