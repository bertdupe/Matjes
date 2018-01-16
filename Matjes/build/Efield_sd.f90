      module m_energyfield
      use m_derived_types
      interface Efield_sd
       module procedure Efield_sd_serial
      end interface Efield_sd

      contains

      subroutine Efield_sd_serial(j,spin,shape_spin,tableNN,masque,indexNN,h_ext,my_lattice)
      use m_energy
      use m_parameters, only : i_DM,i_biq,i_four,i_dip,EA
      implicit none
! input
      type(lattice), intent(in) :: my_lattice
      real(kind=8), intent(in) :: spin(:,:,:,:,:),h_ext(:)
      integer, intent(in) :: tableNN(:,:,:,:,:,:),masque(:,:,:,:),indexNN(:,:)
      integer, intent(in) :: shape_spin(:)
      integer, intent(in) :: j
! internals
      integer :: i_x,i_y,i_z,i_m,k,i
      character(len=30) :: fname,toto
      real(kind=8) :: E_int,E_DM,E_xch,E_ani,E_z,E_4,E_biq,E_dip,mu_B

      E_DM=0.0d0
      E_xch=0.0d0
      E_ani=0.0d0
      E_z=0.0d0
      E_int=0.0d0
      E_4=0.0d0
      E_biq=0.0d0
      E_dip=0.0d0

      write(fname,'(I12)') j
      toto=trim(adjustl(fname))
      write(fname,'(a,18a,a)')'Efield_',(toto(i:i),i=1, &
         len_trim(toto)),'.dat'
      OPEN(70,FILE=fname,action='write',status='unknown',form='formatted')

      Write(70,'(a)') '1:Posx   2:Posy   3:Posz   4:E_tot   5:E_xch   6:E_DM   7:E_ani   8:E_z   9:E_4   10:E_biq   11:E_dip'

      do i_m=1,shape_spin(5)
       do i_z=1,shape_spin(4)
        do i_y=1,shape_spin(3)
         do i_x=1,shape_spin(2)

         mu_B=spin(7,i_x,i_y,i_z,i_m)

          E_xch=Exchange(i_x,i_y,i_z,i_m,spin,shape_spin,tableNN,masque,indexNN)
          E_z=Zeeman(i_x,i_y,i_z,i_m,spin,shape_spin,masque,h_ext,mu_B)
          E_ani=anisotropy(i_x,i_y,i_z,i_m,EA,spin,shape_spin)
          if (i_DM) E_DM=DMenergy(i_x,i_y,i_z,i_m,spin,shape_spin,tableNN,masque,indexNN)
          if (i_four) E_4=fourspin(i_x,i_y,i_z,i_m,spin,shape_spin,masque,my_lattice%boundary)
          if (i_biq) E_biq=biquadratic(i_x,i_y,i_z,i_m,spin,shape_spin,tableNN,masque,indexNN)
#ifdef CPP_BRUTDIP
          if (i_dip) E_dip=dipole(i_x,i_y,i_z,i_m,spin,shape_spin,masque,my_lattice%boundary)
#else
          if (i_dip) E_dip=fftdip(i_x,i_y,i_z,i_m)
#endif

          E_int=E_xch+E_DM+E_ani+E_z+E_4+E_biq+E_dip

          Write(70,'(3(f10.4,2x),8(E20.12E3,2x))') (Spin(k,i_x,i_y,i_z,i_m),k=1,3),E_int,E_xch,E_DM,E_ani,E_z,E_4,E_biq,E_dip

         enddo
        enddo
       enddo
      enddo

      close(70)

      end subroutine Efield_sd_serial

      end module
