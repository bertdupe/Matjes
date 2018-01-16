      module m_energy
      contains

      real(kind=8) function Exchange(i_x,i_y,i_z,i_m,spin,shape_spin,tableNN,masque,indexNN)
      use m_efield, only : me,Efield_Jij
      use m_parameters, only : J_ij,J_il,c_ji,J_z,param_N_Nneigh,param_N_Nneigh_il
#ifdef CPP_MPI
      use m_mpi_prop, only : start
#endif
      implicit none
! input
      integer, intent(in) :: shape_spin(:)
      real(kind=8), intent(in) :: spin(:,:,:,:,:)
      integer, intent(in) :: tableNN(:,:,:,:,:,:),masque(:,:,:,:),indexNN(:,:)
! external variable
      real (kind=8) :: E_int
      integer , intent(in) :: i_x,i_y,i_z,i_m
! position of the neighbors along x,y,z and motif
      integer :: v_x,v_y,v_z,v_m
! position of the spins with the ghost
      integer :: ig_x,ig_y,ig_z
! corrdinate of the spin
      integer :: X,Y,Z
! internal variable
      integer :: i,j,avant,k,lay
#ifndef CPP_MPI
      integer, dimension(3) :: start=0
#endif

      E_int=0.0d0
      Exchange=0.0d0
      X=shape_spin(1)-3
      Y=shape_spin(1)-2
      Z=shape_spin(1)-1
      ig_x=i_x-start(1)
      ig_y=i_y-start(2)
      ig_z=i_z-start(3)

      avant=0
! the choosen spin is the spin i_s. It has nn nearest neighbours. The numbers of nearest
! neighbours are stored in n(:). for example n(1)=4 means 4 nearest neighbours
!!! first neighb'
       if (size(J_ij,2).ne.1) then
        lay=i_m
       else
        lay=1
       endif

       k=1
       ! param_N_Nneigh is the highest order of J, that is unequal to zero in "inp"c
       do i=1,param_N_Nneigh
        do j=1,indexNN(i,k)
        v_x=tableNN(1,avant+j,ig_x,ig_y,ig_z,i_m)
        v_y=tableNN(2,avant+j,ig_x,ig_y,ig_z,i_m)
        v_z=tableNN(3,avant+j,ig_x,ig_y,ig_z,i_m)
        v_m=tableNN(4,avant+j,ig_x,ig_y,ig_z,i_m)

        E_int=E_int+(Spin(X,i_x,i_y,i_z,i_m)*Spin(X,v_x,v_y,v_z,v_m)+ &
         Spin(Y,i_x,i_y,i_z,i_m)*Spin(Y,v_x,v_y,v_z,v_m)+ &
         Spin(Z,i_x,i_y,i_z,i_m)*Spin(Z,v_x,v_y,v_z,v_m))* &
         dble(masque(avant+j+1,i_x,i_y,i_z))* &
         (J_ij(i,lay)+me(i)*Efield_Jij(i_x,i_y,i_z))

        enddo
        avant=avant+indexNN(i,k)
       enddo

       if ((shape_spin(5).ne.1).or.(dabs(J_il(1)).gt.1.0d-8)) then
       avant=sum(indexNN(:,1))
       k=2
       do i=1,param_N_Nneigh_il
        do j=1,indexNN(i,k)
        v_x=tableNN(1,avant+j,ig_x,ig_y,ig_z,i_m)
        v_y=tableNN(2,avant+j,ig_x,ig_y,ig_z,i_m)
        v_z=tableNN(3,avant+j,ig_x,ig_y,ig_z,i_m)
        v_m=tableNN(4,avant+j,ig_x,ig_y,ig_z,i_m)

        E_int=E_int+(Spin(X,i_x,i_y,i_z,i_m)*Spin(X,v_x,v_y,v_z,v_m)+ &
         Spin(Y,i_x,i_y,i_z,i_m)*Spin(Y,v_x,v_y,v_z,v_m)+ &
         Spin(Z,i_x,i_y,i_z,i_m)*Spin(Z,v_x,v_y,v_z,v_m))* &
         dble(masque(avant+j+1,i_x,i_y,i_z))*J_il(i)

        enddo
        avant=avant+indexNN(i,k)
       enddo
       endif

       if ((shape_spin(4).ne.1).or.(dabs(J_z(1)).gt.1.0d-8)) then
       avant=sum(indexNN(:,1))+sum(indexNN(:,2))
       k=3
       do i=1,count(dabs(J_z)>1.0d-8)
        do j=1,indexNN(i,k)
        v_x=tableNN(1,avant+j,ig_x,ig_y,ig_z,i_m)
        v_y=tableNN(2,avant+j,ig_x,ig_y,ig_z,i_m)
        v_z=tableNN(3,avant+j,ig_x,ig_y,ig_z,i_m)
        v_m=tableNN(4,avant+j,ig_x,ig_y,ig_z,i_m)

        E_int=E_int+(Spin(X,i_x,i_y,i_z,i_m)*Spin(X,v_x,v_y,v_z,v_m)+ &
         Spin(Y,i_x,i_y,i_z,i_m)*Spin(Y,v_x,v_y,v_z,v_m)+ &
         Spin(Z,i_x,i_y,i_z,i_m)*Spin(Z,v_x,v_y,v_z,v_m))* &
         dble(masque(avant+j+1,i_x,i_y,i_z))*J_z(i)

        enddo
        avant=avant+indexNN(i,k)
       enddo
       endif

      Exchange=E_int*c_Ji
      end function Exchange

!Zeeman energy
      real(kind=8) function Zeeman(i_x,i_y,i_z,i_m,spin,shape_spin,masque,H_ext,mu_B)
      implicit none
! input
      integer, intent(in) :: shape_spin(:)
      real(kind=8), intent(in) :: spin(:,:,:,:,:)
      real(kind=8), intent(in) ::H_ext(:),mu_B
      integer, intent(in) :: masque(:,:,:,:)
! external variable
      integer , intent(in) :: i_x,i_y,i_z,i_m
! components of the spin
      integer :: X,Y,Z,M

      X=shape_spin(1)-3
      Y=shape_spin(1)-2
      Z=shape_spin(1)-1
      M=shape_spin(1)

      Zeeman=-mu_B*(H_ext(1)*Spin(X,i_x,i_y,i_z,i_m)+H_ext(2)*Spin(Y,i_x,i_y,i_z,i_m)+&
      H_ext(3)*Spin(Z,i_x,i_y,i_z,i_m))*Spin(M,i_x,i_y,i_z,i_m)*dble(masque(1,i_x,i_y,i_z))

      end function Zeeman

!DM energy
      real(kind=8) function DMenergy(i_x,i_y,i_z,i_m,spin,shape_spin,tableNN,masque,indexNN)
      use m_vector, only: TripleProduct,cross
      use m_parameters, only : DM,DM_vector,c_DM
#ifdef CPP_MPI
      use m_mpi_prop, only : start
#endif
      implicit none
! input
      real(kind=8), intent(in) :: spin(:,:,:,:,:)
      integer, intent(in) :: tableNN(:,:,:,:,:,:),masque(:,:,:,:),indexNN(:,:),shape_spin(:)
! external variable
      real(kind=8) :: E_int
      integer , intent(in) :: i_x,i_y,i_z,i_m
! position of the neighbors along x,y,z and motif
      integer :: v_x,v_y,v_z,v_m,ig_x,ig_y,ig_z
! coordinates of the components of the spins
      integer :: X,Y,Z,M
!internal variable
      integer :: avant,i,j
#ifndef CPP_MPI
      integer, dimension(3) :: start=0
#endif

      avant=0

      E_int=0.0d0
      DMenergy=0.0d0
      X=shape_spin(1)-3
      Y=shape_spin(1)-2
      Z=shape_spin(1)-1
      M=shape_spin(1)

      ig_x=i_x-start(1)
      ig_y=i_y-start(2)
      ig_z=i_z-start(3)

      if ((size(DM,2).eq.1).and.(i_m.eq.2)) then
       DMenergy=0.0d0
       return
      endif

      do i=1,count(dabs(DM(:,i_m))>1.0d-8)
       do j=1,indexNN(i,1)
        v_x=tableNN(1,avant+j,ig_x,ig_y,ig_z,i_m)
        v_y=tableNN(2,avant+j,ig_x,ig_y,ig_z,i_m)
        v_z=tableNN(3,avant+j,ig_x,ig_y,ig_z,i_m)
        v_m=tableNN(4,avant+j,ig_x,ig_y,ig_z,i_m)

        E_int=E_int+DM(i,i_m)*TripleProduct(DM_vector(avant+j,:,i_m), &
         Spin(X:Z,i_x,i_y,i_z,i_m),Spin(X:Z,v_x,v_y,v_z,v_m))* &
         dble(masque(avant+j+1,i_x,i_y,i_z))
       enddo
       avant=avant+indexNN(i,1)
      enddo

      DMenergy=c_DM*E_int

      end function DMenergy
!!! Anisotropy
      real(kind=8) function anisotropy(i_x,i_y,i_z,i_m,axis,spin,shape_spin)
      use m_vector, only : norm
      use m_parameters, only : D_ani,c_ani,EA
      implicit none
! input
      real(kind=8), intent(in) :: spin(:,:,:,:,:)
      integer, intent(in) :: shape_spin(:)
! external variable
      real(kind=8) :: E_int
      real(kind=8) , intent(in) :: axis(3)
      integer , intent(in) :: i_x,i_y,i_z,i_m
! components of the spin
      integer :: X,Y,Z,M
! internal variable
      integer :: i
      real(kind=8) :: dumy(3)


      E_int=0.0d0
      anisotropy=0.0d0
      X=shape_spin(1)-3
      Y=shape_spin(1)-2
      Z=shape_spin(1)-1
      M=shape_spin(1)

      dumy=axis/norm(EA)
      do i=1,count(dabs(D_ani)>1.0d-8)
      E_int=D_ani(i)*(dumy(1)*spin(X,i_x,i_y,i_z,i_m)+dumy(2)*spin(Y,i_x,i_y,i_z,i_m)+ &
      dumy(3)*spin(Z,i_x,i_y,i_z,i_m))**2+E_int
      enddo

      anisotropy=c_ani*E_int
      end function anisotropy

! biquadratic energy
      real(kind=8) function biquadratic(i_x,i_y,i_z,i_m,spin,shape_spin,tableNN,masque,indexNN)
      use m_parameters, only : c_JB,J_B
#ifdef CPP_MPI
      use m_mpi_prop, only : start
#endif
      implicit none
! input
      integer, intent(in) :: shape_spin(:)
      real(kind=8), intent(in) :: spin(:,:,:,:,:)
      integer, intent(in) :: tableNN(:,:,:,:,:,:),masque(:,:,:,:),indexNN(:,:)
! external variable
      real(kind=8) :: E_int
      integer , intent(in) :: i_x,i_y,i_z,i_m
! internal variable
      integer :: j
! position of the neighbors along x,y,z and motif
      integer :: v_x,v_y,v_z,v_m,ig_x,ig_y,ig_z
! components of the magnetic moments
      integer :: X,Y,Z,M
#ifndef CPP_MPI
      integer, dimension(3) :: start=0
#endif
      E_int=0.0d0
      biquadratic=0.0d0
      X=shape_spin(1)-3
      Y=shape_spin(1)-2
      Z=shape_spin(1)-1
      M=shape_spin(1)

      ig_x=i_x-start(1)
      ig_y=i_y-start(2)
      ig_z=i_z-start(3)

        do j=1,indexNN(1,1)
        v_x=tableNN(1,j,ig_x,ig_y,ig_z,i_m)
        v_y=tableNN(2,j,ig_x,ig_y,ig_z,i_m)
        v_z=tableNN(3,j,ig_x,ig_y,ig_z,i_m)
        v_m=tableNN(4,j,ig_x,ig_y,ig_z,i_m)

       E_int=E_int+(Spin(X,i_x,i_y,i_z,i_m)*Spin(X,v_x,v_y,v_z,v_m)+ &
        Spin(Y,i_x,i_y,i_z,i_m)*Spin(Y,v_x,v_y,v_z,v_m)+ &
        Spin(Z,i_x,i_y,i_z,i_m)*Spin(Z,v_x,v_y,v_z,v_m))**2*dble(masque(j+1,i_x,i_y,i_z))
        enddo

      biquadratic=c_JB*J_B*E_int*dble(masque(1,i_x,i_y,i_z))
      end function biquadratic
!4 spin term
      real(kind=8) function fourspin(i_x,i_y,i_z,i_m,spin,shape_spin,masque,Periodic_log)
      use m_sym_utils, only : corners
      use m_parameters, only :c_Ki,K_1
      implicit none
! input
      real(kind=8), intent(in) :: spin(:,:,:,:,:)
      integer, intent(in) :: masque(:,:,:,:),shape_spin(:)
      logical, intent(in) :: Periodic_log(:)
! external variable
      real(kind=8) :: E_int,E_local
      integer , intent(in) :: i_x,i_y,i_z,i_m
! components of the spin
      integer :: X,Y,Z,M
! internal variable
      integer :: k,ipu(3),ipv(3),ipuv(3)

      E_int=0.0d0
      fourspin=0.0d0
      X=shape_spin(1)-3
      Y=shape_spin(1)-2
      Z=shape_spin(1)-1
      M=shape_spin(1)

      k=1
      do while (k.lt.size(corners,1))
       E_local=0.0d0

! take care of the non-periodic boundary condition
       if (.not.Periodic_log(1)) then
        if (((i_x-1+corners(k+2,1)).lt.0).or.((i_x-1+corners(k+2,1)).ge.shape_spin(2))) then
         k=k+3
         cycle
        endif

        if (((i_x-1+corners(k,1)).lt.0).or.((i_x-1+corners(k,1)).ge.shape_spin(2))) then
         k=k+3
         cycle
        endif

        if (((i_x-1+corners(k+1,1)).lt.0).or.((i_x-1+corners(k+1,1)).ge.shape_spin(2))) then
         k=k+3
         cycle
        endif
       endif

       if (.not.Periodic_log(2)) then
        if (((i_y-1+corners(k+2,2)).lt.0).or.((i_y-1+corners(k+2,2)).ge.shape_spin(3))) then
         k=k+3
         cycle
        endif

        if (((i_y-1+corners(k,2)).lt.0).or.((i_y-1+corners(k,2)).ge.shape_spin(3))) then
         k=k+3
         cycle
        endif

        if (((i_y-1+corners(k+1,2)).lt.0).or.((i_y-1+corners(k+1,2)).ge.shape_spin(3))) then
         k=k+3
         cycle
        endif
       endif

       ipuv=(/mod(i_x-1+corners(k+2,1)+shape_spin(2),shape_spin(2))+1, &
        mod(i_y+corners(k+2,2)+shape_spin(3)-1,shape_spin(3))+1, &
        i_z/)

       ipu=(/mod(i_x-1+corners(k,1)+shape_spin(2),shape_spin(2))+1, &
        mod(i_y+corners(k,2)+shape_spin(3)-1,shape_spin(3))+1, &
        i_z/)

       ipv=(/mod(i_x-1+corners(k+1,1)+shape_spin(2),shape_spin(2))+1, &
        mod(i_y+corners(k+1,2)+shape_spin(3)-1,shape_spin(3))+1, &
        i_z/)

       E_local=E_local+(Spin(X,i_x,i_y,i_z,i_m)*Spin(X,ipu(1),ipu(2),ipu(3),i_m)+ &
         Spin(Y,i_x,i_y,i_z,i_m)*Spin(Y,ipu(1),ipu(2),ipu(3),i_m)+ &
         Spin(Z,i_x,i_y,i_z,i_m)*Spin(Z,ipu(1),ipu(2),ipu(3),i_m))*&
         (Spin(X,ipv(1),ipv(2),ipv(3),i_m)*Spin(X,ipuv(1),ipuv(2),ipuv(3),i_m)+Spin(Y,ipv(1),ipv(2),ipv(3),i_m)* &
         Spin(Y,ipuv(1),ipuv(2),ipuv(3),i_m)+Spin(Z,ipv(1),ipv(2),ipv(3),i_m)*Spin(Z,ipuv(1),ipuv(2),ipuv(3),i_m))

       E_local=E_local+(Spin(X,i_x,i_y,i_z,i_m)*Spin(X,ipv(1),ipv(2),ipv(3),i_m)+ &
         Spin(Y,i_x,i_y,i_z,i_m)*Spin(Y,ipv(1),ipv(2),ipv(3),i_m)+ &
         Spin(Z,i_x,i_y,i_z,i_m)*Spin(Z,ipv(1),ipv(2),ipv(3),i_m))*&
         (Spin(X,ipu(1),ipu(2),ipu(3),i_m)*Spin(X,ipuv(1),ipuv(2),ipuv(3),i_m)+Spin(Y,ipu(1),ipu(2),ipu(3),i_m)* &
         Spin(Y,ipuv(1),ipuv(2),ipuv(3),i_m)+Spin(Z,ipu(1),ipu(2),ipu(3),i_m)*Spin(Z,ipuv(1),ipuv(2),ipuv(3),i_m))

       E_local=E_local-(Spin(X,i_x,i_y,i_z,i_m)*Spin(X,ipuv(1),ipuv(2),ipuv(3),i_m)+Spin(Y,i_x,i_y,i_z,i_m)* &
         Spin(Y,ipuv(1),ipuv(2),ipuv(3),i_m)+Spin(Z,i_x,i_y,i_z,i_m)*Spin(Z,ipuv(1),ipuv(2),ipuv(3),i_m))* &
         (Spin(X,ipu(1),ipu(2),ipu(3),i_m)*Spin(X,ipv(1),ipv(2),ipv(3),i_m)+Spin(Y,ipu(1),ipu(2),ipu(3),i_m)* &
         Spin(Y,ipv(1),ipv(2),ipv(3),i_m)+Spin(Z,ipu(1),ipu(2),ipu(3),i_m)*Spin(Z,ipv(1),ipv(2),ipv(3),i_m))

       E_local=E_local*dble(masque(1,i_x,i_y,i_z)*masque(1,ipu(1),ipu(2),ipu(3))* &
       masque(1,ipv(1),ipv(2),ipv(3))*masque(1,ipuv(1),ipuv(2),ipuv(3)))
       k=k+3
       E_int=E_int+E_local
      enddo

      fourspin=c_Ki*K_1*E_int

      end function fourspin

#ifdef CPP_BRUTDIP

! Dipole Dipole interaction
      real(kind=8) function dipole(i_x,i_y,i_z,i_m,spin,shape_spin,masque,Periodic_log)
      use m_constants, only : pi
      use m_vector, only : norm
      implicit none
! input
      real(kind=8), intent(in) :: spin(:,:,:,:,:)
      integer, intent(in) :: masque(:,:,:,:),shape_spin(:)
      logical, intent(in) :: Periodic_log(:)
! external variable
      integer , intent(in) :: i_x,i_y,i_z,i_m
! coordinate of the magnetic moments
      integer :: Rx,Ry,Rz
      integer :: X,Y,Z,M
! internal variable
      integer :: j_x,j_y,j_z,j_m,nmag
      real(kind=8) :: rc(3),step,ss
      real(kind=8), parameter :: alpha=6.74582d-7
! alpha=mu0*muB/a**3; a in nm. The result is an energy in eV

      dipole=0.0d0
      nmag=shape_spin(5)
      Rx=shape_spin(1)-6
      Ry=shape_spin(1)-5
      Rz=shape_spin(1)-4
      X=shape_spin(1)-3
      Y=shape_spin(1)-2
      Z=shape_spin(1)-1
      M=shape_spin(1)

! the choosen spin is the spin i_s. It has nn nearest neighbours. The numbers of nearest
! neighbours are stored in n(:). for example n(1)=4 means 4 nearest neighbours
!!! first neighb

       do j_m=1,nmag
        do j_z=1,shape_spin(4)
         do j_y=1,shape_spin(3)
          do j_x=1,shape_spin(2)

          if (masque(1,j_x,j_y,j_z).eq.0) cycle

          if (all(Periodic_log)) then
            rc=spin(Rx:Rz,j_x,j_y,j_z,j_m)-spin(Rx:Rz,i_x,i_y,i_z,i_m)
          else
            stop "not implemented"
          endif

          ss=norm(rc)
          if (ss.lt.1.0d-3) cycle
          rc=rc/ss

          step=dot_product(rc,spin(X:Z,i_x,i_y,i_z,i_m))*dot_product(rc,spin(X:Z,j_x,j_y,j_z,j_m))

          dipole=dipole+(dot_product(spin(X:Z,j_x,j_y,j_z,j_m),spin(X:Z,i_x,i_y,i_z,i_m))-3.0d0*step)/ss**3 &
       &   *spin(M,i_x,i_y,i_z,i_m)*spin(M,j_x,j_y,j_z,j_m)

          enddo
         enddo
        enddo
       enddo

      dipole=dipole/pi(4.0d0)*0.5d0*alpha

      end function dipole
#endif

!stoner energy
      real(kind=8) function stoner(i_x,i_y,i_z,i_m,spin,masque,tableNN,indexNN)
      use m_parameters, only : Ist
      implicit none
! input
      real(kind=8), intent(in) :: spin(:,:,:,:,:)
      integer, intent(in) :: tableNN(:,:,:,:,:,:),masque(:,:,:,:),indexNN(:,:)
! external variable
      integer , intent(in) :: i_x,i_y,i_z,i_m
!internals
      integer :: v_x,v_y,v_z,v_m
      integer :: j,avant
      real(kind=8) :: E_int

      E_int=0.0d0
      stoner=0.0d0

      avant=0

        do j=1,indexNN(1,1)
        v_x=tableNN(1,avant+j,i_x,i_y,i_z,i_m)
        v_y=tableNN(2,avant+j,i_x,i_y,i_z,i_m)
        v_z=tableNN(3,avant+j,i_x,i_y,i_z,i_m)
        v_m=tableNN(4,avant+j,i_x,i_y,i_z,i_m)

        E_int=E_int+Spin(7,i_x,i_y,i_z,i_m)*(Spin(4,i_x,i_y,i_z,i_m)*Spin(4,v_x,v_y,v_z,v_m)+ &
         Spin(5,i_x,i_y,i_z,i_m)*Spin(5,v_x,v_y,v_z,v_m)+ &
         Spin(6,i_x,i_y,i_z,i_m)*Spin(6,v_x,v_y,v_z,v_m))**2* &
         dble(masque(avant+j+1,i_x,i_y,i_z))

        enddo

      stoner=E_int*Ist*dble(masque(1,i_x,i_y,i_z))/2.0d0

      end function stoner

#ifndef CPP_BRUTDIP
!dipole dipole energy from fft
      real(kind=8) function fftdip(i_x,i_y,i_z,i_m)
      use m_setup_dipole, only : mmatrix,mcomplex,hreal,hcomplex,Nfftx,Nffty,Nfftz,ntensor &
     & ,rtrans,ctrans,planrtoc,planctor
      use m_lattice, only : spin
      use m_rw_lattice, only : dim_lat
      use m_constants, only : mu_B
      use m_fft
      implicit none
      integer, intent(in) :: i_x,i_y,i_z,i_m
      !dummy
      integer :: i,j,k,l

      ! update the local changed spin
      mmatrix(:,i_x,i_y,i_z)=spin(4:6,i_x,i_y,i_z,1)*spin(7,i_x,i_y,i_z,1)

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

#ifdef CPP_DEBUG
      do k=1,dim_lat(3)
       do l=1,3
       write(*,*) 'coordinate ', l, 'layer ', k
        do j=1,dim_lat(2)
       write(*,*) (hreal(l,i,j,k),i=1,dim_lat(1))
        enddo
       enddo
       pause
      enddo
#endif

      fftdip=(Spin(4,i_x,i_y,i_z,1)*hreal(1,i_x,i_y,i_z)+Spin(5,i_x,i_y,i_z,1)*hreal(2,i_x,i_y,i_z)+ &
       Spin(6,i_x,i_y,i_z,1)*hreal(3,i_x,i_y,i_z))*mu_B

      end function fftdip
#endif

      end module
