      module m_local_energy
      interface local_energy
       module procedure local_energy_decompose
       module procedure local_energy_av
       module procedure local_energy_MC
      end interface local_energy
      contains

      function local_energy_MC(i_DM,i_four,i_biq,i_dip,i_stone,EA,Ilat, &
          & spin,shape_spin,tableNN,shape_tableNN,masque,shape_masque,indexNN,shape_index,h_ext)
      use m_energy
      implicit none
! input
      integer, intent(in) :: shape_index(2),shape_spin(5),shape_tableNN(6),shape_masque(4)
      real(kind=8), intent(in) :: spin(shape_spin(1),shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5)),h_ext(3),EA(3)
      integer, intent(in) :: tableNN(shape_tableNN(1),shape_tableNN(2),shape_tableNN(3),shape_tableNN(4),shape_tableNN(5),shape_tableNN(6))
      integer, intent(in) :: masque(shape_masque(1),shape_masque(2),shape_masque(3),shape_masque(4)),indexNN(shape_index(1),shape_index(2))
      integer, intent(in) :: Ilat(:)
      logical, intent(in) :: i_DM,i_four,i_biq,i_dip,i_stone
! ouput
      real(kind=8) :: local_energy_MC(8)
! internal
      real(kind=8) :: E_int(8),mu_B
      integer :: i_x,i_y,i_z,i_m

      E_int=0.0d0
      local_energy_MC=0.0d0

      i_x=Ilat(1)
      i_y=Ilat(2)
      i_z=Ilat(3)
      i_m=Ilat(4)

      mu_B=spin(7,i_x,i_y,i_z,i_m)

      E_int(1)=2.0d0*Exchange(i_x,i_y,i_z,i_m,spin,shape_spin,tableNN,masque,indexNN)
      E_int(2)=Zeeman(i_x,i_y,i_z,i_m,spin,shape_spin,masque,h_ext,mu_B)
      E_int(3)=anisotropy(i_x,i_y,i_z,i_m,EA,spin,shape_spin,masque)
      if (i_DM) E_int(5)=2.0d0*DMenergy(i_x,i_y,i_z,i_m,spin,shape_spin,tableNN,masque,indexNN)
      if (i_four) E_int(4)=4.0d0*fourspin(i_x,i_y,i_z,i_m,spin,shape_spin,masque)
      if (i_biq) E_int(6)=2.0d0*biquadratic(i_x,i_y,i_z,i_m,spin,shape_spin,tableNN,masque,indexNN)
#ifdef CPP_BRUTDIP
      if (i_dip) E_int(7)=dipole(i_x,i_y,i_z,i_m,spin,shape_spin,masque)
#else
      if (i_dip) E_int(7)=fftdip(i_x,i_y,i_z,i_m)
#endif
      if (i_stone) E_int(8)=stoner(i_x,i_y,i_z,i_m,spin,masque,tableNN,indexNN)

      local_energy_MC=E_int

      end function local_energy_MC

      function local_energy_decompose(E_dim,i_DM,i_four,i_biq,i_dip,EA,i_x,i_y,i_z,i_m, &
          & spin,shape_spin,tableNN,shape_tableNN,masque,shape_masque,indexNN,shape_index,h_ext)
      use m_energy
      implicit none
! input
      integer, intent(in) :: shape_index(2),shape_spin(5),shape_tableNN(6),shape_masque(4)
      real(kind=8), intent(in) :: spin(shape_spin(1),shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5)),h_ext(3),EA(3)
      integer, intent(in) :: tableNN(shape_tableNN(1),shape_tableNN(2),shape_tableNN(3),shape_tableNN(4),shape_tableNN(5),shape_tableNN(6))
      integer, intent(in) :: masque(shape_masque(1),shape_masque(2),shape_masque(3),shape_masque(4)),indexNN(shape_index(1),shape_index(2))
      integer, intent(in) :: i_x,i_y,i_z,i_m
      logical, intent(in) :: i_DM,i_four,i_biq,i_dip
      real(kind=8), intent(in) :: E_dim(:)
! ouput
      real(kind=8) :: local_energy_decompose(8)
! internal
      real(kind=8) :: E_int(8),mu_B

      E_int=0.0d0
      mu_B=spin(7,i_x,i_y,i_z,i_m)

      E_int(1)=Exchange(i_x,i_y,i_z,i_m,spin,shape_spin,tableNN,masque,indexNN)
      E_int(2)=Zeeman(i_x,i_y,i_z,i_m,spin,shape_spin,masque,h_ext,mu_B)
      E_int(3)=anisotropy(i_x,i_y,i_z,i_m,EA,spin,shape_spin,masque)
      if (i_DM) E_int(5)=DMenergy(i_x,i_y,i_z,i_m,spin,shape_spin,tableNN,masque,indexNN)
      if (i_four) E_int(4)=fourspin(i_x,i_y,i_z,i_m,spin,shape_spin,masque)
      if (i_biq) E_int(6)=biquadratic(i_x,i_y,i_z,i_m,spin,shape_spin,tableNN,masque,indexNN)
#ifdef CPP_BRUTDIP
      if (i_dip) E_int(7)=dipole(i_x,i_y,i_z,i_m,spin,shape_spin,masque)
#else
      if (i_dip) E_int(7)=fftdip(i_x,i_y,i_z,i_m)
#endif

      local_energy_decompose=E_int

      end function local_energy_decompose

      real(kind=8) function local_energy_av(E_dim,i_DM,i_four,i_biq,i_dip,EA,i_x,i_y,i_z,i_m, &
          & spin,shape_spin,tableNN,shape_tableNN,masque,shape_masque,indexNN,shape_index,h_ext)
      use m_energy
      implicit none
! input
      integer, intent(in) :: shape_index(2),shape_spin(5),shape_tableNN(6),shape_masque(4)
      real(kind=8), intent(in) :: spin(shape_spin(1),shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5)),h_ext(3),EA(3)
      integer, intent(in) :: tableNN(shape_tableNN(1),shape_tableNN(2),shape_tableNN(3),shape_tableNN(4),shape_tableNN(5),shape_tableNN(6))
      integer, intent(in) :: masque(shape_masque(1),shape_masque(2),shape_masque(3),shape_masque(4)),indexNN(shape_index(1),shape_index(2))
      integer, intent(in) :: i_x,i_y,i_z,i_m
      logical, intent(in) :: i_DM,i_four,i_biq,i_dip
      real(kind=8), intent(in) :: E_dim
! internal
      real(kind=8) :: E_int,mu_B

      E_int=0.0d0
      mu_B=spin(7,i_x,i_y,i_z,i_m)

      E_int=Exchange(i_x,i_y,i_z,i_m,spin,shape_spin,tableNN,masque,indexNN)
      E_int=E_int+Zeeman(i_x,i_y,i_z,i_m,spin,shape_spin,masque,h_ext,mu_B)
      E_int=E_int+anisotropy(i_x,i_y,i_z,i_m,EA,spin,shape_spin,masque)
      if (i_DM) E_int=E_int+DMenergy(i_x,i_y,i_z,i_m,spin,shape_spin,tableNN,masque,indexNN)
      if (i_four) E_int=E_int+fourspin(i_x,i_y,i_z,i_m,spin,shape_spin,masque)
      if (i_biq) E_int=E_int+biquadratic(i_x,i_y,i_z,i_m,spin,shape_spin,tableNN,masque,indexNN)
#ifdef CPP_BRUTDIP
      if (i_dip) E_int=E_int+dipole(i_x,i_y,i_z,i_m,spin,shape_spin,masque)
#else
      if (i_dip) E_int=E_int+fftdip(i_x,i_y,i_z,i_m)
#endif

      local_energy_av=E_int

      end function local_energy_av

      end module m_local_energy