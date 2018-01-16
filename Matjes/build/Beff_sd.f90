      module m_eval_Beff
      interface calculate_Beff
       module procedure normal
      end interface calculate_Beff
      contains
! subroutine that calculates the field
! dE/DM
!
!
!
!--------------------------------------------------------------
! for normal
      subroutine normal(i_DM,i_four,i_biq,i_dip,EA,i_x,i_y,i_z,i_m,B, &
         &   spin,shape_spin,indexNN,shape_index,masque,shape_masque,tableNN,shape_tableNN,h_ext, &
         my_lattice)
      use m_fieldeff
      use m_derived_types
      implicit none
! input variable
      logical, intent(in) :: i_DM,i_four,i_biq,i_dip
      integer, intent(in) :: i_x,i_y,i_z,i_m
      integer, intent(in) :: shape_index(2),shape_spin(5),shape_tableNN(6),shape_masque(4)
      integer, intent(in) :: tableNN(shape_tableNN(1),shape_tableNN(2),shape_tableNN(3),shape_tableNN(4),shape_tableNN(5),shape_tableNN(6))
      integer, intent(in) :: masque(shape_masque(1),shape_masque(2),shape_masque(3),shape_masque(4))
      integer, intent(in) :: indexNN(shape_index(1),shape_index(2))
      real(kind=8), intent(in) :: spin(shape_spin(1),shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5)),h_ext(3),EA(3)
      type(lattice), intent(in) :: my_lattice
! output of the function
      real(kind=8), intent(out) :: B(:)
! internals

       if (masque(1,i_x,i_y,i_z).eq.0) then
         B=0.0d0
         return
       endif

        B=Bexch(i_x,i_y,i_z,i_m,spin,shape_spin,indexNN,shape_index,masque,shape_masque,tableNN,shape_tableNN) &
     &   +BZ(i_x,i_y,i_z,i_m,spin,shape_spin,h_ext)+Bani(i_x,i_y,i_z,i_m,EA,spin,shape_spin) &
     &   +efield(i_x,i_y,i_z,i_m,spin,shape_spin,indexNN,shape_index,masque,shape_masque,tableNN,shape_tableNN)
        if (i_DM) B=B+Bdm(i_x,i_y,i_z,i_m,spin,shape_spin,indexNN,shape_index,masque,shape_masque,tableNN,shape_tableNN)
        if (i_four) B=B+Bfour(i_x,i_y,i_z,i_m,spin,shape_spin,masque,shape_masque,my_lattice)
        if (i_biq) B=B+Bbiqd(i_x,i_y,i_z,i_m,spin,shape_spin,indexNN,shape_index,masque,shape_masque,tableNN,shape_tableNN)
        if (i_dip) B=B+Bdip(i_x,i_y,i_z,i_m,spin,shape_spin)

      end subroutine normal

      end module m_eval_Beff
