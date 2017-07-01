      module m_relaxtyp
      contains
! functions that relaxes the spins with respect of the dE/dM
! in one case, the spins are choosen in the direction of -dE/DM so energy diminishes
! in the second case, the spins are choosen in the direction of +dE/DM so energy increases
!
      function underrelax(spin_in,i_x,i_y,i_z,i_m,spin,shape_spin,indexNN,shape_index,masque,shape_masque,tableNN,shape_tableNN,h_ext)
      use m_parameters, only : EA,i_biq,i_four,i_DM
      use m_vector, only : cross
      use m_fieldeff
      implicit none
! external variable
      real(kind=8), dimension(3), intent(in) :: spin_in,h_ext
      integer, intent(in) :: i_x,i_y,i_z,i_m
      integer, intent(in) :: shape_index(2),shape_spin(5),shape_tableNN(6),shape_masque(4)
      integer, intent(in) :: tableNN(shape_tableNN(1),shape_tableNN(2),shape_tableNN(3),shape_tableNN(4),shape_tableNN(5),shape_tableNN(6))
      integer, intent(in) :: masque(shape_masque(1),shape_masque(2),shape_masque(3),shape_masque(4))
      integer, intent(in) :: indexNN(shape_index(1),shape_index(2))
      real(kind=8), intent(in) :: spin(shape_spin(1),shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5))
! value of the function
      real(kind=8), dimension(3) :: underrelax
!internal variable
      real(kind=8), dimension(3) ::S_int
      real(kind=8) :: norm,dumy(3)

      S_int=0.0d0
      S_int=S_int+Bexch(i_x,i_y,i_z,i_m,spin,shape_spin,indexNN,shape_index,masque,shape_masque,tableNN,shape_tableNN)
      S_int=S_int+BZ(i_x,i_y,i_z,i_m,spin,shape_spin,h_ext)
      if (i_DM) S_int=S_int+Bdm(i_x,i_y,i_z,i_m,spin,shape_spin,indexNN,shape_index,masque,shape_masque,tableNN,shape_tableNN)
      S_int=S_int+Bani(i_x,i_y,i_z,i_m,EA,spin,shape_spin,masque,shape_masque)
      S_int=S_int+efield(i_x,i_y,i_z,i_m,spin,shape_spin,indexNN,shape_index,masque,shape_masque,tableNN,shape_tableNN)
      if (i_biq) S_int=S_int+Bbiqd(i_x,i_y,i_z,i_m,spin,shape_spin,indexNN,shape_index,masque,shape_masque,tableNN,shape_tableNN)
      if (i_four) S_int=S_int+Bfour(i_x,i_y,i_z,i_m,spin,shape_spin,masque,shape_masque)

      norm=sqrt(S_int(1)**2+S_int(2)**2+S_int(3)**2)
! Calculation of the new spin
!      norm=dsqrt((S_int(1)+spin_in(1))**2+(S_int(2)+spin_in(2))**2+(S_int(3)+spin_in(3))**2)
!      underrelax=(S_int+spin_in)/norm

      if (norm.gt.1.0d-8) then
       dumy=cross(spin_in,S_int)
       S_int=spin_in-cross(spin_in,dumy)
       norm=sqrt(S_int(1)**2+S_int(2)**2+S_int(3)**2)
       underrelax=S_int/norm
      else
       underrelax=spin_in
      endif

      end function underrelax

      function overrelax(spin_in,i_x,i_y,i_z,i_m,spin,shape_spin,indexNN,shape_index,masque,shape_masque,tableNN,shape_tableNN,h_ext)
      use m_parameters, only : EA,i_biq,i_four,i_DM
      use m_vector, only : cross
      use m_fieldeff
      implicit none
! external variable
      real(kind=8), dimension(3), intent(in) :: spin_in,h_ext
      integer, intent(in) :: i_x,i_y,i_z,i_m
      integer, intent(in) :: shape_index(2),shape_spin(5),shape_tableNN(6),shape_masque(4)
      integer, intent(in) :: tableNN(shape_tableNN(1),shape_tableNN(2),shape_tableNN(3),shape_tableNN(4),shape_tableNN(5),shape_tableNN(6))
      integer, intent(in) :: masque(shape_masque(1),shape_masque(2),shape_masque(3),shape_masque(4))
      integer, intent(in) :: indexNN(shape_index(1),shape_index(2))
      real(kind=8), intent(in) :: spin(shape_spin(1),shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5))
! value of the function
      real(kind=8), dimension(3) :: overrelax
!internal variable
      real(kind=8), dimension(3) ::S_int
      real(kind=8) :: norm,dumy(3)

      S_int=0.0d0

      S_int=S_int+Bexch(i_x,i_y,i_z,i_m,spin,shape_spin,indexNN,shape_index,masque,shape_masque,tableNN,shape_tableNN)
      S_int=S_int+BZ(i_x,i_y,i_z,i_m,spin,shape_spin,h_ext)
      if (i_DM) S_int=S_int+Bdm(i_x,i_y,i_z,i_m,spin,shape_spin,indexNN,shape_index,masque,shape_masque,tableNN,shape_tableNN)
      S_int=S_int+Bani(i_x,i_y,i_z,i_m,EA,spin,shape_spin,masque,shape_masque)
      S_int=S_int+efield(i_x,i_y,i_z,i_m,spin,shape_spin,indexNN,shape_index,masque,shape_masque,tableNN,shape_tableNN)
      if (i_biq) S_int=S_int+Bbiqd(i_x,i_y,i_z,i_m,spin,shape_spin,indexNN,shape_index,masque,shape_masque,tableNN,shape_tableNN)
      if (i_four) S_int=S_int+Bfour(i_x,i_y,i_z,i_m,spin,shape_spin,masque,shape_masque)

      norm=sqrt(S_int(1)**2+S_int(2)**2+S_int(3)**2)
! Calculation of the new spin
!      norm=dsqrt((-S_int(1)+spin_in(1))**2+(-S_int(2)+spin_in(2))**2+(-S_int(3)+spin_in(3))**2)
!      overrelax=(-S_int+spin_in)/norm

      if (norm.gt.1.0d-8) then
       dumy=cross(spin_in,S_int)
       S_int=spin_in-cross(spin_in,dumy)
       norm=sqrt(S_int(1)**2+S_int(2)**2+S_int(3)**2)
       overrelax=-S_int/norm
      else
       overrelax=-spin_in
      endif

      end function overrelax
      end module
