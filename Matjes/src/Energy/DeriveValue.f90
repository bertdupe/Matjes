! initiate the value of the total energy for the starting configuration
! E is the decomposition of the total energy and is an out.
      subroutine DeriveValue(spin,shape_spin,tableNN,shape_tableNN,masque,shape_masque,indexNN,shape_index,E, &
                 &    i_four,i_dm,i_biq,i_dip,i_stone,EA,h_ext,my_lattice)
      use m_derived_types
      use m_total_energy
      Implicit none
      type(lattice), intent(in) :: my_lattice
      integer, intent(in) :: shape_index(2),shape_spin(5),shape_tableNN(6),shape_masque(4)
      real(kind=8), intent(in) :: spin(shape_spin(1),shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5))
      integer, intent(in) :: tableNN(shape_tableNN(1),shape_tableNN(2),shape_tableNN(3),shape_tableNN(4),shape_tableNN(5),shape_tableNN(6))
      integer, intent(in) :: masque(shape_masque(1),shape_masque(2),shape_masque(3),shape_masque(4))
      integer, intent(in) :: indexNN(shape_index(1),shape_index(2))
      real(kind=8), intent(in) :: EA(3),h_ext(3)
      logical, intent(in) :: i_four,i_dm,i_biq,i_dip,i_stone
      real(kind=8), intent(out) :: E(8)
!     slope variable

      E=0.0d0

      E(1)=total_Exchange(spin)
      E(2)=total_Zeeman(spin,h_ext)
      E(3)=total_anisotropy(spin)
!      if (i_four) E(4)=total_fourspin(spin,shape_spin,masque,shape_masque,my_lattice%boundary)
      if (i_dm) E(5)=total_DMenergy(spin)
!      if (i_biq) E(6)=total_biquadratic(spin,shape_spin,tableNN,shape_tableNN,masque,shape_masque,indexNN,shape_index)

!      if (i_stone) E(8)=total_stoner(spin,shape_spin,masque,shape_masque,tableNN,shape_tableNN,indexNN,shape_index)

      write(6,'(a)') 'Reinitialization of the total Energies'
      write(6,'(a,/)') 'Energy of the spin configuration (in eV)'
      write(6,'(a,4x,f16.8)') 'Exchange energy', E(1)
      write(6,'(a,4x,f16.8)') 'Zeeman energy', E(2)
      write(6,'(a,4x,f16.8)') 'Anisotropy energy', E(3)
      if (i_four) write(6,'(a,4x,f16.8)') '4-spin energy', E(4)
      if (i_dm) write(6,'(a,4x,f16.8)') 'DMI energy', E(5)
      if (i_biq) write(6,'(a,4x,f16.8)') 'Biquadratic energy', E(6)
      if (i_dip) write(6,'(a,4x,f16.8)') 'Dipole-dipole energy', E(7)
      write(6,'(a,4x,f16.8)') 'Total energy', sum(E)

      end subroutine DeriveValue
