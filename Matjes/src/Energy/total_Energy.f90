      module m_total_energy
! this module contains the routines that calculate the total energy for all the energy terms
! input
! spin: the spin lattice. The table is not split
! shape_spin: the size of the lattice table
! shape_masque: same for the masque
! shape_...: same for table ...
! tableNN: the table of neighbors for each of the spins. Be carefull, this table is split and the indices of the spins
! i_x, i_y and i_z do not refer to the correct spins. All the coordinates have to be translated.
! masque: gives the presence or absence of the spins on site i_x, i_y, i_z. This table is also split so the coordinates
! have to be translated
! indexNN: gives the number of atoms in the first, second... shells
! EA: gives the easy axis
! all the rest are read from modules so look into the module
      contains

      real(kind=8) function total_Exchange(spin)
!      use m_efield, only : me,Efield_Jij
#ifdef CPP_MPI
      use m_make_box, only : Xstart,Ystart,Zstart
      use m_mpi_prop, only : start
#endif
      implicit none
      real(kind=8), intent(in) :: spin(:,:,:,:,:)
! input
      total_Exchange=0.0d0

      end function total_Exchange


!Zeeman energy
      real(kind=8) function total_Zeeman(spin,h_ext)
#ifdef CPP_MPI
      use m_make_box, only : Xstart,Xstop,Ystart,Ystop,Zstart,Zstop
      use m_mpi_prop, only : start
#endif
      implicit none
! input
      real(kind=8), intent(in) :: spin(:,:,:,:,:),h_ext(3)
! external variable

      total_Zeeman=0.0d0

      total_Zeeman=-sum(spin(1,:,:,:,:)*h_ext(1)+spin(2,:,:,:,:)*h_ext(2)+spin(3,:,:,:,:)*h_ext(3))

      end function total_Zeeman

!DM energy
      real(kind=8) function total_DMenergy(spin)
      use m_vector, only: TripleProduct
#ifdef CPP_MPI
      use m_make_box, only : Xstart,Ystart,Zstart
      use m_mpi_prop, only : start
#endif
      implicit none
! input
      real(kind=8), intent(in) :: spin(:,:,:,:,:)
! external variable


!         E_int=E_int+DM(i,i_m)*TripleProduct(DM_vector(avant+j,:,i_m), &
!          Spin(4:6,i_x,i_y,i_z,i_m),Spin(4:6,v_x,v_y,v_z,v_m))* &
!          dble(masque(avant+j+1,i_x,i_y,i_z)*masque(1,i_x,i_y,i_z))


      total_DMenergy=0.0d0

      end function total_DMenergy
!!! Anisotropy
      real(kind=8) function total_anisotropy(spin)
#ifdef CPP_MPI
      use m_make_box, only : Xstart,Xstop,Ystart,Ystop,Zstart,Zstop
#endif
      implicit none
! input
      real(kind=8), intent(in) :: spin(:,:,:,:,:)
! internal variable

      total_anisotropy=0.0d0

      end function total_anisotropy

! biquadratic energy
!     real(kind=8) function total_biquadratic(spin,shape_spin,tableNN,shape_tableNN,masque,shape_masque,indexNN,shape_index)
!      use m_parameters, only : c_JB,J_B
!#ifdef CPP_MPI
!      use m_make_box, only : Xstart,Ystart,Zstart
!      use m_mpi_prop, only : start
!#endif
!      implicit none
!! input
!      integer, intent(in) :: shape_index(2),shape_spin(5),shape_tableNN(6),shape_masque(4)
!      real(kind=8), intent(in) :: spin(shape_spin(1),shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5))
!      integer, intent(in) :: tableNN(shape_tableNN(1),shape_tableNN(2),shape_tableNN(3),shape_tableNN(4),shape_tableNN(5),shape_tableNN(6))
!      integer, intent(in) :: masque(shape_masque(1),shape_masque(2),shape_masque(3),shape_masque(4))
!      integer, intent(in) :: indexNN(shape_index(1),shape_index(2))
!! value of the function
!      real(kind=8) :: E_int
!      integer :: i_x,i_y,i_z,i_m
!! internal variable
!      integer :: j
!! position of the neighbors along x,y,z and motif
!      integer :: v_x,v_y,v_z,v_m,ig_x,ig_y,ig_z
!#ifndef CPP_MPI
!      integer :: Xstart,Ystart,Zstart
!      integer, dimension(3) :: start=0
!
!      Xstart=1
!      Ystart=1
!      Zstart=1
!#endif
!
!      E_int=0.0d0
!      total_biquadratic=0.0d0
!! the choosen spin is the spin i_s. It has nn nearest neighbours. The numbers of nearest
!! neighbours are stored in n(:). for example n(1)=4 means 4 nearest neighbours
!!!! first neighb
!#ifdef CPP_OPENMP
!!$OMP parallel DO REDUCTION(+:E_int) private(i_x,i_y,i_z,i_m,j,v_x,v_y,v_z,v_m) default(shared)
!#endif
!      do j=1,indexNN(1,1)
!      do i_m=1,shape_tableNN(6)
!       do ig_z=1,shape_tableNN(5)
!        do ig_y=1,shape_tableNN(4)
!         do ig_x=1,shape_tableNN(3)
!
!          i_x=ig_x+start(1)
!          i_y=ig_y+start(2)
!          i_z=ig_z+start(3)
!
!          v_x=tableNN(1,j,ig_x,ig_y,ig_z,i_m)
!          v_y=tableNN(2,j,ig_x,ig_y,ig_z,i_m)
!          v_z=tableNN(3,j,ig_x,ig_y,ig_z,i_m)
!          v_m=tableNN(4,j,ig_x,ig_y,ig_z,i_m)
!
!       E_int=E_int+J_B*(Spin(4,i_x,i_y,i_z,i_m)*Spin(4,v_x,v_y,v_z,v_m)+ &
!        Spin(5,i_x,i_y,i_z,i_m)*Spin(5,v_x,v_y,v_z,v_m)+&
!        Spin(6,i_x,i_y,i_z,i_m)*Spin(6,v_x,v_y,v_z,v_m))**2*dble(masque(j+1,i_x,i_y,i_z)*masque(1,i_x,i_y,i_z))
!         enddo
!        enddo
!       enddo
!      enddo
!      enddo
!#ifdef CPP_OPENMP
!!$OMP end parallel do
!#endif
!      total_biquadratic=c_JB*E_int
!      end function total_biquadratic
!!4 spin term
!      real(kind=8) function total_fourspin(spin,shape_spin,masque,shape_masque,Periodic_log)
!      use m_parameters, only :c_Ki,K_1
!      use m_sym_utils, only : corners
!#ifdef CPP_MPI
!      use m_make_box, only : Xstart,Xstop,Ystart,Ystop,Zstart,Zstop
!#endif
!      implicit none
!! input
!      integer, intent(in) :: shape_spin(5),shape_masque(4)
!      real(kind=8), intent(in) :: spin(shape_spin(1),shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5))
!      integer, intent(in) :: masque(shape_masque(1),shape_masque(2),shape_masque(3),shape_masque(4))
!      logical, intent(in) :: Periodic_log(3)
!! external variable
!      real(kind=8) :: E_int,E_local
!! internal variable
!      integer :: k,ipu(3),ipv(3),ipuv(3),i_x,i_y,i_z,i_m
!#ifndef CPP_MPI
!      integer :: Xstop,Xstart,Ystop,Ystart,Zstop,Zstart
!
!      Xstop=shape_spin(2)
!      Xstart=1
!      Ystop=shape_spin(3)
!      Ystart=1
!      Zstop=shape_spin(4)
!      Zstart=1
!#endif
!
!      E_int=0.0d0
!
!#ifdef CPP_OPENMP
!!$OMP parallel DO REDUCTION(+:E_int) private(i_x,i_y,i_z,i_m,ipuv,ipu,ipv,k) default(shared)
!#endif
!      do i_m=1,shape_spin(5)
!      do i_z=Zstart,Zstop
!       do i_y=Ystart,Ystop
!        do i_x=Xstart,Xstop
!
!       k=1
!       do while (k.lt.size(corners,1))
!        E_local=0.0d0
!
!! take care of the non-periodic boundary condition
!       if (.not.Periodic_log(1)) then
!        if (((i_x-1+corners(k+2,1)).lt.0).or.((i_x-1+corners(k+2,1)).ge.shape_spin(2))) then
!         k=k+3
!         cycle
!        endif
!
!        if (((i_x-1+corners(k,1)).lt.0).or.((i_x-1+corners(k,1)).ge.shape_spin(2))) then
!         k=k+3
!         cycle
!        endif
!
!        if (((i_x-1+corners(k+1,1)).lt.0).or.((i_x-1+corners(k+1,1)).ge.shape_spin(2))) then
!         k=k+3
!         cycle
!        endif
!       endif
!
!       if (.not.Periodic_log(2)) then
!        if (((i_y-1+corners(k+2,2)).lt.0).or.((i_y-1+corners(k+2,2)).ge.shape_spin(3))) then
!         k=k+3
!         cycle
!        endif
!
!        if (((i_y-1+corners(k,2)).lt.0).or.((i_y-1+corners(k,2)).ge.shape_spin(3))) then
!         k=k+3
!         cycle
!        endif
!
!        if (((i_y-1+corners(k+1,2)).lt.0).or.((i_y-1+corners(k+1,2)).ge.shape_spin(3))) then
!         k=k+3
!         cycle
!        endif
!       endif
!
!       ipuv=(/mod(i_x-1+corners(k+2,1)+shape_spin(2),shape_spin(2))+1, &
!        mod(i_y+corners(k+2,2)+shape_spin(3)-1,shape_spin(3))+1, &
!        i_z/)
!
!       ipu=(/mod(i_x-1+corners(k,1)+shape_spin(2),shape_spin(2))+1, &
!        mod(i_y+corners(k,2)+shape_spin(3)-1,shape_spin(3))+1, &
!        i_z/)
!
!       ipv=(/mod(i_x-1+corners(k+1,1)+shape_spin(2),shape_spin(2))+1, &
!        mod(i_y+corners(k+1,2)+shape_spin(3)-1,shape_spin(3))+1, &
!        i_z/)
!
!        E_local=E_local+(Spin(4,i_x,i_y,i_z,i_m)*Spin(4,ipu(1),ipu(2),ipu(3),i_m)+ &
!         Spin(5,i_x,i_y,i_z,i_m)*Spin(5,ipu(1),ipu(2),ipu(3),i_m)+ &
!         Spin(6,i_x,i_y,i_z,i_m)*Spin(6,ipu(1),ipu(2),ipu(3),i_m))* &
!         (Spin(4,ipv(1),ipv(2),ipv(3),i_m)*Spin(4,ipuv(1),ipuv(2),ipuv(3),i_m)+ &
!         Spin(5,ipv(1),ipv(2),ipv(3),i_m)*Spin(5,ipuv(1),ipuv(2),ipuv(3),i_m)+ &
!         Spin(6,ipv(1),ipv(2),ipv(3),i_m)*Spin(6,ipuv(1),ipuv(2),ipuv(3),i_m))
!
!        E_local=E_local+(Spin(4,i_x,i_y,i_z,i_m)*Spin(4,ipv(1),ipv(2),ipv(3),i_m)+ &
!         Spin(5,i_x,i_y,i_z,i_m)*Spin(5,ipv(1),ipv(2),ipv(3),i_m)+ &
!         Spin(6,i_x,i_y,i_z,i_m)*Spin(6,ipv(1),ipv(2),ipv(3),i_m))* &
!         (Spin(4,ipu(1),ipu(2),ipu(3),i_m)*Spin(4,ipuv(1),ipuv(2),ipuv(3),i_m)+ &
!         Spin(5,ipu(1),ipu(2),ipu(3),i_m)*Spin(5,ipuv(1),ipuv(2),ipuv(3),i_m)+ &
!         Spin(6,ipu(1),ipu(2),ipu(3),i_m)*Spin(6,ipuv(1),ipuv(2),ipuv(3),i_m))
!
!        E_local=E_local-(Spin(4,i_x,i_y,i_z,i_m)*Spin(4,ipuv(1),ipuv(2),ipuv(3),i_m)+ &
!         Spin(5,i_x,i_y,i_z,i_m)*Spin(5,ipuv(1),ipuv(2),ipuv(3),i_m)+ &
!         Spin(6,i_x,i_y,i_z,i_m)*Spin(6,ipuv(1),ipuv(2),ipuv(3),i_m))* &
!         (Spin(4,ipu(1),ipu(2),ipu(3),i_m)*Spin(4,ipv(1),ipv(2),ipv(3),i_m)+ &
!         Spin(5,ipu(1),ipu(2),ipu(3),i_m)*Spin(5,ipv(1),ipv(2),ipv(3),i_m)+ &
!         Spin(6,ipu(1),ipu(2),ipu(3),i_m)*Spin(6,ipv(1),ipv(2),ipv(3),i_m))
!
!        E_local=E_local*K_1*dble(masque(1,i_x,i_y,i_z)*masque(1,ipu(1),ipu(2),ipu(3))* &
!         masque(1,ipv(1),ipv(2),ipv(3))*masque(1,ipuv(1),ipuv(2),ipuv(3)))
!        E_int=E_int+E_local
!        k=k+3
!       enddo
!      enddo
!      enddo
!      enddo
!      enddo
!#ifdef CPP_OPENMP
!!$OMP end parallel do
!#endif
!      total_fourspin=c_Ki*E_int
!
!      end function total_fourspin

! Dipole Dipole interaction
      real(kind=8) function total_dipole(spin,shape_spin,Periodic_log)
      use m_constants, only : pi
      use m_vector, only : norm
#ifdef CPP_MPI
      use m_make_box, only : Xstart,Xstop,Ystart,Ystop,Zstart,Zstop
#endif
      implicit none
! input
      integer, intent(in) :: shape_spin(5)
      real(kind=8), intent(in) :: spin(shape_spin(1),shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5))
      logical, intent(in) :: Periodic_log(3)
! external variable
      integer :: i_x,i_y,i_z,j_x,j_y,j_z,i_m,j_m,nmag
! internal variable
      real(kind=8) :: rc(3),ss
      real(kind=8), parameter :: alpha=6.74582d-7
#ifndef CPP_MPI
      integer :: Xstop,Xstart,Ystop,Ystart,Zstop,Zstart

      Xstop=shape_spin(2)
      Xstart=1
      Ystop=shape_spin(3)
      Ystart=1
      Zstop=shape_spin(4)
      Zstart=1
#endif

      total_dipole=0.0d0
      nmag=shape_spin(5)
! the choosen spin is the spin i_s. It has nn nearest neighbours. The numbers of nearest
! neighbours are stored in n(:). for example n(1)=4 means 4 nearest neighbours
!!! first neighb

#ifdef CPP_OPENMP
!$OMP parallel do reduction(+:total_dipole) private(i_x,i_y,i_z,j_x,j_y,j_z,j,rc,ss,i_m,j_m) default(shared)
#endif

       do i_m=1,nmag
        do i_z=Zstart,Zstop
         do i_y=Ystart,Ystop
          do i_x=Xstart,Xstop

         if (spin(7,i_x,i_y,i_z,i_m).eq.0) cycle

       do j_m=1,nmag
        do j_z=1,shape_spin(4)
         do j_y=1,shape_spin(3)
          do j_x=1,shape_spin(2)

          if (all(periodic_log)) then
            rc=spin(1:3,j_x,j_y,j_z,j_m)-spin(1:3,i_x,i_y,i_z,i_m)
          endif

         ss=norm(rc)
         if (ss.lt.1.0d-3) cycle
         rc=rc/ss

         total_dipole=total_dipole+(dot_product(spin(4:6,j_x,j_y,j_z,j_m),spin(4:6,i_x,i_y,i_z,i_m))-3.0d0* &
          dot_product(rc,spin(4:6,i_x,i_y,i_z,i_m))*dot_product(rc,spin(4:6,j_x,j_y,j_z,j_m))/ss**3)* &
          spin(7,j_x,j_y,j_z,j_m)

          enddo
         enddo
        enddo
       enddo

          enddo
         enddo
        enddo
       enddo
#ifdef CPP_OPENMP
!$OMP end parallel do
#endif
      total_dipole=total_dipole/pi(4.0d0)*0.5d0*alpha

      write(*,*) total_dipole

      end function total_dipole

! total stoner energy
!      real(kind=8) function total_stoner(spin,shape_spin,masque,shape_masque,tableNN,shape_tableNN,indexNN,shape_index)
!      use m_parameters, only : Ist
!      implicit none
!! input
!      integer, intent(in) :: shape_index(2),shape_spin(5),shape_tableNN(6),shape_masque(4)
!      real(kind=8), intent(in) :: spin(shape_spin(1),shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5))
!      integer, intent(in) :: tableNN(shape_tableNN(1),shape_tableNN(2),shape_tableNN(3),shape_tableNN(4),shape_tableNN(5),shape_tableNN(6))
!      integer, intent(in) :: masque(shape_masque(1),shape_masque(2),shape_masque(3),shape_masque(4))
!      integer, intent(in) :: indexNN(shape_index(1),shape_index(2))
!! external variable
!      integer :: i_x,i_y,i_z,i_m
!!internals
!      integer :: v_x,v_y,v_z,v_m
!      integer :: j,avant
!      real(kind=8) :: E_int
!
!      E_int=0.0d0
!      avant=0
!
!      do i_m=1,shape_spin(5)
!      do i_z=1,shape_spin(4)
!       do i_y=1,shape_spin(3)
!        do i_x=1,shape_spin(2)
!
!        do j=1,indexNN(1,1)
!        v_x=tableNN(1,avant+j,i_x,i_y,i_z,i_m)
!        v_y=tableNN(2,avant+j,i_x,i_y,i_z,i_m)
!        v_z=tableNN(3,avant+j,i_x,i_y,i_z,i_m)
!        v_m=tableNN(4,avant+j,i_x,i_y,i_z,i_m)
!
!        E_int=E_int+Spin(7,i_x,i_y,i_z,i_m)*(Spin(4,i_x,i_y,i_z,i_m)*Spin(4,v_x,v_y,v_z,v_m)+ &
!         Spin(5,i_x,i_y,i_z,i_m)*Spin(5,v_x,v_y,v_z,v_m)+ &
!         Spin(6,i_x,i_y,i_z,i_m)*Spin(6,v_x,v_y,v_z,v_m))**2* &
!         dble(masque(avant+j+1,i_x,i_y,i_z))
!        enddo
!
!        enddo
!       enddo
!      enddo
!      enddo
!
!      total_stoner=E_int*Ist*dble(masque(1,i_x,i_y,i_z))/2.0d0
!
!      end function total_stoner

#ifndef CPP_BRUTDIP
! total energy of the fft dipole
      real(kind=8) function total_fftdip()
      use m_setup_dipole, only : mmatrix,mcomplex,hreal,hcomplex,Nfftx,Nffty,Nfftz,ntensor &
     & ,rtrans,ctrans,planrtoc,planctor
      use m_lattice, only : spin
      use m_rw_lattice, only : dim_lat
      use m_constants, only : mu_B
      use m_fft
      implicit none
      !dummy
      integer :: i,j,k,l

      total_fftdip=0.0d0

      ! FFT with libraries depending on your choice

      call fft(Nfftx,Nffty,Nfftz,mmatrix,mcomplex,rtrans,ctrans,planrtoc)

      ! calculate field and depolarising energy

      do k=1,Nfftz !z
       do j=1,Nffty !y
        do i=1,Nfftx !x

      hcomplex(1,i,j,k)=ntensor(1,i,j,k)*mcomplex(1,i,j,k)+ntensor(4,i,j,k)*mcomplex(2,i,j,k)+ &
       ntensor(6,i,j,k)*mcomplex(3,i,j,k)
      hcomplex(2,i,j,k)=ntensor(4,i,j,k)*mcomplex(1,i,j,k)+ntensor(2,i,j,k)*mcomplex(2,i,j,k)+ &
       ntensor(5,i,j,k)*mcomplex(3,i,j,k)
      hcomplex(3,i,j,k)=ntensor(6,i,j,k)*mcomplex(1,i,j,k)+ntensor(5,i,j,k)*mcomplex(2,i,j,k)+ &
       ntensor(3,i,j,k)*mcomplex(3,i,j,k)

        enddo
       enddo
      enddo

      call fft(Nfftx,Nffty,Nfftz,hcomplex,hreal,-1,rtrans,ctrans,planctor)

      do k=1,dim_lat(3)
       do j=1,dim_lat(2)
        do i=1,dim_lat(1)
         total_fftdip=total_fftdip+(Spin(4,i,j,k,1)*hreal(1,i,j,k)+ &
         Spin(5,i,j,k,1)*hreal(2,i,j,k)+ &
         Spin(6,i,j,k,1)*hreal(3,i,j,k))
        enddo
       enddo
      enddo

      total_fftdip=total_fftdip*mu_B/2.0d0

      end function total_fftdip
#endif
      end module
