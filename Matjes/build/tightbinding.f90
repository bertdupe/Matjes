      subroutine tightbinding(spin,shape_spin)
      implicit none
!      integer, intent(in) :: shape_spin(5)
!      real(kind=8), intent(in) :: spin(shape_spin(1),shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5))

      integer, intent(in) :: shape_spin
      real(kind=8), intent(in) :: spin(shape_spin)
      real(kind=8) :: toto(shape_spin)
      write(6,'(a)') 'your now in the tight-binding code'

      toto=spin

      end subroutine tightbinding
