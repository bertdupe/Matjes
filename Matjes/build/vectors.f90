!!!!  Simple functions that can deal with vector vector calculations
!!!!
!!!!
      #include "math.h"
      module m_vectors
      interface norm
       module procedure norm_real
       module procedure norm_int
      end interface norm
      contains


!
!     Cross product
!
      function cross_real(a, b)
      implicit none
      real(Kind=DEF_DBL_PREC), dimension(3) :: cross_real
      real(Kind=DEF_DBL_PREC), intent(in) :: a(:), b(:)

      cross_real(1) = a(2) * b(3) - a(3) * b(2)
      cross_real(2) = a(3) * b(1) - a(1) * b(3)
      cross_real(3) = a(1) * b(2) - a(2) * b(1)
      end function cross_real

      function cross_int(a, b)
      implicit none
      real(Kind=DEF_DBL_PREC), dimension(3) :: cross_int
      integer, intent(in) :: a(:), b(:)

      cross_int(1) = dble(a(2) * b(3) - a(3) * b(2))
      cross_int(2) = dble(a(3) * b(1) - a(1) * b(3))
      cross_int(3) = dble(a(1) * b(2) - a(2) * b(1))
      end function cross_int

      end module m_vectors
