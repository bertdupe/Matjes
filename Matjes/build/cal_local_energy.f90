module m_local_energy
use m_derived_types

      interface local_energy_MC
       module procedure loc_en_MC
      end interface local_energy_MC

      interface local_energy
       module procedure local_energy_decompose
       module procedure local_energy_av
       module procedure local_energy_pointer
      end interface local_energy

private
public :: local_energy_MC,local_energy,local_energy_pointer

contains

      function loc_en_MC(i_dip,Ilat, &
          & spin,shape_spin,tableNN,shape_tableNN,masque,shape_masque,indexNN,shape_index,h_ext, &
          & my_lattice,Hamiltonian)
      use m_energy
      implicit none
! input
      integer, intent(in) :: shape_index(2),shape_spin(5),shape_tableNN(6),shape_masque(4)
      real(kind=8), intent(in) :: spin(:,:,:,:,:),h_ext(3)
      integer, intent(in) :: tableNN(shape_tableNN(1),shape_tableNN(2),shape_tableNN(3),shape_tableNN(4),shape_tableNN(5),shape_tableNN(6))
      integer, intent(in) :: masque(shape_masque(1),shape_masque(2),shape_masque(3),shape_masque(4)),indexNN(shape_index(1),shape_index(2))
      integer, intent(in) :: Ilat(:)
      logical, intent(in) :: i_dip
      type(lattice), intent(in) :: my_lattice
      type(Coeff_Ham), intent(in) :: Hamiltonian
! ouput
      real(kind=8) :: loc_en_MC(8)
! internal
      real(kind=8) :: E_int(8),mu_B
      integer :: i_x,i_y,i_z,i_m
      logical :: Periodic_log(3)

      E_int=0.0d0
      loc_en_MC=0.0d0

      i_x=Ilat(1)
      i_y=Ilat(2)
      i_z=Ilat(3)
      i_m=Ilat(4)

      mu_B=sqrt(spin(1,i_x,i_y,i_z,i_m)**2+spin(2,i_x,i_y,i_z,i_m)**2+spin(3,i_x,i_y,i_z,i_m)**2)
      Periodic_log=my_lattice%boundary

      if (Hamiltonian%i_exch) E_int(1)=2.0d0*Exchange(i_x,i_y,i_z,i_m,spin,shape_spin,tableNN,masque,indexNN,Hamiltonian%c_Ji,Hamiltonian%exchange(:,:,1))
      E_int(2)=Zeeman(i_x,i_y,i_z,i_m,spin,shape_spin,masque,h_ext,mu_B)
      if (Hamiltonian%i_ani) E_int(3)=anisotropy(i_x,i_y,i_z,i_m,spin,shape_spin,Hamiltonian%ani,Hamiltonian%c_ani)
      if (Hamiltonian%i_DM) E_int(5)=2.0d0*DMenergy(i_x,i_y,i_z,i_m,spin,shape_spin,tableNN,masque,indexNN,Hamiltonian%c_DM,Hamiltonian%DMI)
      if (Hamiltonian%i_four) E_int(4)=4.0d0*fourspin(i_x,i_y,i_z,i_m,spin,shape_spin,masque,Periodic_log,Hamiltonian%c_Ki,Hamiltonian%fours)
      if (Hamiltonian%i_biq) E_int(6)=2.0d0*biquadratic(i_x,i_y,i_z,i_m,spin,shape_spin,tableNN,masque,indexNN,Hamiltonian%c_Jb,Hamiltonian%biq)
#ifdef CPP_BRUTDIP
      if (i_dip) E_int(7)=dipole(i_x,i_y,i_z,i_m,spin,shape_spin,masque,Periodic_log)
#else
      if (i_dip) E_int(7)=fftdip(i_x,i_y,i_z,i_m)
#endif
!      if (i_stone) E_int(8)=stoner(i_x,i_y,i_z,i_m,spin,masque,tableNN,indexNN)

      loc_en_MC=E_int

      end function loc_en_MC

      subroutine local_energy_decompose(E_int,i_dip,i_x,i_y,i_z,i_m, &
          & spin,shape_spin,tableNN,shape_tableNN,masque,shape_masque,indexNN,shape_index,h_ext, &
          & my_lattice,Hamiltonian)
      use m_energy
      implicit none
! input
      integer, intent(in) :: shape_index(2),shape_spin(5),shape_tableNN(6),shape_masque(4)
      real(kind=8), intent(in) :: spin(shape_spin(1),shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5)),h_ext(3)
      integer, intent(in) :: tableNN(shape_tableNN(1),shape_tableNN(2),shape_tableNN(3),shape_tableNN(4),shape_tableNN(5),shape_tableNN(6))
      integer, intent(in) :: masque(shape_masque(1),shape_masque(2),shape_masque(3),shape_masque(4)),indexNN(shape_index(1),shape_index(2))
      integer, intent(in) :: i_x,i_y,i_z,i_m
      logical, intent(in) :: i_dip
      type(lattice), intent(in) :: my_lattice
      type(Coeff_Ham), intent(in) :: Hamiltonian
! ouput
      real(kind=8), intent(out) :: E_int(:)
! internal
      real(kind=8) :: mu_B
      logical :: Periodic_log(3)

      E_int=0.0d0
      mu_B=sqrt(spin(1,i_x,i_y,i_z,i_m)**2+spin(2,i_x,i_y,i_z,i_m)**2+spin(3,i_x,i_y,i_z,i_m)**2)
      Periodic_log=my_lattice%boundary

      if (masque(1,i_x,i_y,i_z).eq.0) return

      if (Hamiltonian%i_exch) E_int(1)=Exchange(i_x,i_y,i_z,i_m,spin,shape_spin,tableNN,masque,indexNN,Hamiltonian%c_Ji,Hamiltonian%exchange(:,:,1))
      E_int(2)=Zeeman(i_x,i_y,i_z,i_m,spin,shape_spin,masque,h_ext,mu_B)
      if (Hamiltonian%i_ani) E_int(3)=anisotropy(i_x,i_y,i_z,i_m,spin,shape_spin,Hamiltonian%ani,Hamiltonian%c_ani)
      if (Hamiltonian%i_DM) E_int(5)=DMenergy(i_x,i_y,i_z,i_m,spin,shape_spin,tableNN,masque,indexNN,Hamiltonian%c_DM,Hamiltonian%DMI)
      if (Hamiltonian%i_four) E_int(4)=fourspin(i_x,i_y,i_z,i_m,spin,shape_spin,masque,Periodic_log,Hamiltonian%c_Ki,Hamiltonian%fours)
      if (Hamiltonian%i_biq) E_int(6)=biquadratic(i_x,i_y,i_z,i_m,spin,shape_spin,tableNN,masque,indexNN,Hamiltonian%c_Jb,Hamiltonian%biq)
#ifdef CPP_BRUTDIP
      if (i_dip) E_int(7)=dipole(i_x,i_y,i_z,i_m,spin,shape_spin,masque,Periodic_log)
#else
      if (i_dip) E_int(7)=fftdip(i_x,i_y,i_z,i_m)
#endif

      end subroutine local_energy_decompose

      subroutine local_energy_av(E_int,i_dip,i_x,i_y,i_z,i_m, &
          & spin,shape_spin,tableNN,shape_tableNN,masque,shape_masque,indexNN,shape_index,h_ext, &
          & my_lattice,Hamiltonian)
      use m_energy
      implicit none
! input
      integer, intent(in) :: shape_index(2),shape_spin(5),shape_tableNN(6),shape_masque(4)
      real(kind=8), intent(in) :: spin(shape_spin(1),shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5)),h_ext(3)
      integer, intent(in) :: tableNN(shape_tableNN(1),shape_tableNN(2),shape_tableNN(3),shape_tableNN(4),shape_tableNN(5),shape_tableNN(6))
      integer, intent(in) :: masque(shape_masque(1),shape_masque(2),shape_masque(3),shape_masque(4)),indexNN(shape_index(1),shape_index(2))
      integer, intent(in) :: i_x,i_y,i_z,i_m
      logical, intent(in) :: i_dip
      type(lattice), intent(in) :: my_lattice
      type(Coeff_Ham), intent(in) :: Hamiltonian
! ouput
      real(kind=8), intent(out) :: E_int
! internal
      real(kind=8) :: mu_B
      logical :: Periodic_log(3)

      E_int=0.0d0
      mu_B=sqrt(spin(1,i_x,i_y,i_z,i_m)**2+spin(2,i_x,i_y,i_z,i_m)**2+spin(3,i_x,i_y,i_z,i_m)**2)
      Periodic_log=my_lattice%boundary

      if (masque(1,i_x,i_y,i_z).eq.0) return

      if (Hamiltonian%i_exch) E_int=Exchange(i_x,i_y,i_z,i_m,spin,shape_spin,tableNN,masque,indexNN,Hamiltonian%c_Ji,Hamiltonian%exchange(:,:,1))
      E_int=E_int+Zeeman(i_x,i_y,i_z,i_m,spin,shape_spin,masque,h_ext,mu_B)
      if (Hamiltonian%i_ani) E_int=E_int+anisotropy(i_x,i_y,i_z,i_m,spin,shape_spin,Hamiltonian%ani,Hamiltonian%c_ani)
      if (Hamiltonian%i_DM) E_int=E_int+DMenergy(i_x,i_y,i_z,i_m,spin,shape_spin,tableNN,masque,indexNN,Hamiltonian%c_DM,Hamiltonian%DMI)
      if (Hamiltonian%i_four) E_int=E_int+fourspin(i_x,i_y,i_z,i_m,spin,shape_spin,masque,Periodic_log,Hamiltonian%c_Ki,Hamiltonian%fours)
      if (Hamiltonian%i_biq) E_int=E_int+biquadratic(i_x,i_y,i_z,i_m,spin,shape_spin,tableNN,masque,indexNN,Hamiltonian%c_Jb,Hamiltonian%biq)
#ifdef CPP_BRUTDIP
      if (i_dip) E_int=E_int+dipole(i_x,i_y,i_z,i_m,spin,shape_spin,masque,Periodic_log)
#else
      if (i_dip) E_int=E_int+fftdip(i_x,i_y,i_z,i_m)
#endif

      end subroutine local_energy_av

subroutine local_energy_pointer(E_int,iomp,spin,E_line)
use m_energy_commons
use m_external_fields, only : ext_field
implicit none
! input
type(point_shell_mode), intent(in) :: spin
type(point_shell_Operator), intent(in) :: E_line
integer, intent(in) :: iomp
! ouput
real(kind=8), intent(out) :: E_int
! internal
integer :: i,N

N=size(spin%shell)
E_int=0.0d0

do i=1,N

   E_int=E_int+dot_product( spin%shell(1)%w , matmul(E_line%shell(i)%Op_loc,spin%shell(i)%w) )

enddo

!#ifdef CPP_BRUTDIP
!      if (i_dip) E_int=E_int+dipole(i_x,i_y,i_z,i_m,spin,shape_spin,masque,Periodic_log)
!#else
!      if (i_dip) E_int=E_int+fftdip(i_x,i_y,i_z,i_m)
!#endif

end subroutine local_energy_pointer

end module m_local_energy
