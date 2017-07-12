      module m_energyfield
      interface Efield_sd
       module procedure Efield_sd_serial
       module procedure Efield_sd_para
      end interface Efield_sd

      contains

      subroutine Efield_sd_serial(j,spin,shape_spin,tableNN,masque,indexNN,h_ext)
      use m_energy
      use m_parameters, only : i_DM,i_biq,i_four,i_dip,EA
      implicit none
! input
      real(kind=8), intent(in) :: spin(:,:,:,:,:),h_ext(:)
      integer, intent(in) :: tableNN(:,:,:,:,:,:),masque(:,:,:,:),indexNN(:,:)
      integer, intent(in) :: shape_spin(:)
      integer, intent(in) :: j
! internals
      integer :: iomp,i_x,i_y,i_z,i_m,k,i
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
          E_ani=anisotropy(i_x,i_y,i_z,i_m,EA,spin,shape_spin,masque)
          if (i_DM) E_DM=DMenergy(i_x,i_y,i_z,i_m,spin,shape_spin,tableNN,masque,indexNN)
          if (i_four) E_4=fourspin(i_x,i_y,i_z,i_m,spin,shape_spin,masque)
          if (i_biq) E_biq=biquadratic(i_x,i_y,i_z,i_m,spin,shape_spin,tableNN,masque,indexNN)
#ifdef CPP_BRUTDIP
          if (i_dip) E_dip=dipole(i_x,i_y,i_z,i_m,spin,shape_spin,masque)
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

! writting the different energy contribution of the energetics with mpi_write commands

      subroutine Efield_sd_para(j,spin,shape_spin,tableNN,shape_tableNN,masque,indexNN,h_ext,irank,start,isize,mpi_comm)
      use m_energy
      use m_parameters, only : i_DM,i_biq,i_four,i_dip,EA
#ifdef CPP_MPI
      use m_mpi
#endif
      implicit none
! input
      real(kind=8), intent(in) :: spin(:,:,:,:,:),h_ext(:)
      integer, intent(in) :: tableNN(:,:,:,:,:,:),masque(:,:,:,:),indexNN(:,:),shape_spin(:),shape_tableNN(:)
      integer, intent(in) :: j,isize,irank,start(:),mpi_comm
! internals
      integer :: iomp,i_x,i_y,i_z,i_m,k,i,i_gx,i_gy,i_gz
      character(len=30) :: fname,toto
      real(kind=8) :: E_density(11,shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5)),mu_B
! mpi variable

      E_density=0.0d0

      do i_m=1,shape_tableNN(6)
       do i_gz=1,shape_tableNN(5)
        do i_gy=1,shape_tableNN(4)
         do i_gx=1,shape_tableNN(3)

          mu_B=spin(7,i_x,i_y,i_z,i_m)

          i_x=i_gx+start(1)
          i_y=i_gy+start(2)
          i_z=i_gz+start(3)

          E_density(1:3,i_x,i_y,i_z,i_m)=spin(1:3,i_x,i_y,i_z,i_m)

          E_density(5,i_x,i_y,i_z,i_m)=Exchange(i_x,i_y,i_z,i_m,spin,shape_spin,tableNN,masque,indexNN)
          E_density(8,i_x,i_y,i_z,i_m)=Zeeman(i_x,i_y,i_z,i_m,spin,shape_spin,masque,h_ext,mu_B)
          E_density(7,i_x,i_y,i_z,i_m)=anisotropy(i_x,i_y,i_z,i_m,EA,spin,shape_spin,masque)
          if (i_DM) E_density(6,i_x,i_y,i_z,i_m)=DMenergy(i_x,i_y,i_z,i_m,spin,shape_spin,tableNN,masque,indexNN)
          if (i_four) E_density(9,i_x,i_y,i_z,i_m)=fourspin(i_x,i_y,i_z,i_m,spin,shape_spin,masque)
          if (i_biq) E_density(10,i_x,i_y,i_z,i_m)=biquadratic(i_x,i_y,i_z,i_m,spin,shape_spin,tableNN,masque,indexNN)
#ifdef CPP_BRUTDIP
          if (i_dip) E_density(11,i_x,i_y,i_z,i_m)=dipole(i_x,i_y,i_z,i_m,spin,shape_spin,masque)
#else
          if (i_dip) E_density(11,i_x,i_y,i_z,i_m)=fftdip(i_x,i_y,i_z,i_m)
#endif

          E_density(4,i_x,i_y,i_z,i_m)=sum(E_density(5:11,i_x,i_y,i_z,i_m))

         enddo
        enddo
       enddo
      enddo

#ifdef CPP_MPI
      E_density=reduce(E_density,11,shape_spin(1),shape_spin(2),shape_spin(3),shape_spin(4),MPI_COMM)
#endif

      if (irank.eq.0) then
       write(fname,'(I12)') j
       toto=trim(adjustl(fname))
       write(fname,'(a,18a,a)')'Efield_',(toto(i:i),i=1,len_trim(toto)),'.dat'
       OPEN(70,FILE=fname,action='write',status='unknown',form='formatted')

       write(70,'(a)') '1:Posx   2:Posy   3:Posz   4:E_tot   5:E_xch   6:E_DM   7:E_ani   8:E_z   9:E_4   10:E_biq   11:E_dip'

       do i_m=1,shape_spin(5)
        do i_z=1,shape_spin(4)
         do i_y=1,shape_spin(3)
          do i_x=1,shape_spin(2)
          Write(70,'(3(f10.4,2x),8(E20.12E3,2x))') (E_density(k,i_x,i_y,i_z,i_m),k=1,11)
          enddo
         enddo
        enddo
       enddo
       close(70)
      endif

      end subroutine Efield_sd_para

      end module
