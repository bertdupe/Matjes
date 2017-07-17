! initiate the value of the total energy for the starting configuration
! E is the decomposition of the total energy and is an out.
      subroutine DeriveValue(N_cell,spin,shape_spin,tableNN,shape_tableNN,masque,shape_masque,indexNN,shape_index,E, &
                 &    i_four,i_dm,i_biq,i_dip,i_stone,EA,h_ext)
      use m_total_energy
      Implicit none
      integer, intent(in) :: N_cell,shape_index(2),shape_spin(5),shape_tableNN(6),shape_masque(4)
      real(kind=8), intent(in) :: spin(shape_spin(1),shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5))
      integer, intent(in) :: tableNN(shape_tableNN(1),shape_tableNN(2),shape_tableNN(3),shape_tableNN(4),shape_tableNN(5),shape_tableNN(6))
      integer, intent(in) :: masque(shape_masque(1),shape_masque(2),shape_masque(3),shape_masque(4))
      integer, intent(in) :: indexNN(shape_index(1),shape_index(2))
      real(kind=8), intent(in) :: EA(3),h_ext(3)
      logical, intent(in) :: i_four,i_dm,i_biq,i_dip,i_stone
      real(kind=8), intent(out) :: E(8)
!     slope variable

      E=0.0d0

      E(1)=total_Exchange(spin,shape_spin,tableNN,shape_tableNN,masque,shape_masque,indexNN,shape_index)
      E(2)=total_Zeeman(spin,shape_spin,masque,shape_masque,h_ext)
      E(3)=total_anisotropy(EA,spin,shape_spin,masque,shape_masque)
      if (i_four) E(4)=total_fourspin(spin,shape_spin,masque,shape_masque)
      if (i_dm) E(5)=total_DMenergy(spin,shape_spin,tableNN,shape_tableNN,masque,shape_masque,indexNN,shape_index)
      if (i_biq) E(6)=total_biquadratic(spin,shape_spin,tableNN,shape_tableNN,masque,shape_masque,indexNN,shape_index)
#ifdef CPP_BRUTDIP
      write(*,*) i_dip
      if (i_dip) E(7)=total_dipole(spin,shape_spin)
#else
      if (i_dip) E(7)=total_fftdip()
#endif
      if (i_stone) E(8)=total_stoner(spin,shape_spin,masque,shape_masque,tableNN,shape_tableNN,indexNN,shape_index)

      end subroutine DeriveValue
