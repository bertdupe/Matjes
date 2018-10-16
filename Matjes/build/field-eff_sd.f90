      module m_fieldeff
      contains
! subfunctions that calculates dE/DM for all the interactions
! input: spin structure with table of neighbors....
!        position of the spin of interest
! output: a 3 vector
!
! May the force be with you
!
!!!! exchange field
      function Bexch(i_x,i_y,i_z,i_m,spin,shape_spin,indexNN,shape_index,masque,shape_masque,tableNN,shape_tableNN,c_ji)
      use m_parameters, only : J_ij,J_il,J_z,param_N_Nneigh,param_N_Nneigh_il
      implicit none
! external variable
      integer , intent(in) :: i_x,i_y,i_z,i_m
      integer, intent(in) :: shape_index(2),shape_spin(5),shape_tableNN(6),shape_masque(4)
      integer, intent(in) :: tableNN(shape_tableNN(1),shape_tableNN(2),shape_tableNN(3),shape_tableNN(4),shape_tableNN(5),shape_tableNN(6))
      integer, intent(in) :: masque(shape_masque(1),shape_masque(2),shape_masque(3),shape_masque(4))
      integer, intent(in) :: indexNN(shape_index(1),shape_index(2))
      real(kind=8), intent(in) :: spin(shape_spin(1),shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5))
      real(kind=8), intent(in) :: c_ji
! value of the function variable
      real(kind=8), dimension(3) :: Bexch
! internal variable
      integer :: i,j,avant,k
! coordinate of the spin
      integer :: X,Y,Z
      real(kind=8), dimension(3) :: B_int
! position of the neighbors along x,y,z and motif
      integer :: v_x,v_y,v_z,v_m,lay

      X=shape_spin(1)-2
      Y=shape_spin(1)-1
      Z=shape_spin(1)-0

      k=1
      B_int=0.0d0
      avant=0

! the choosen spin is the spin i_s. It has nn nearest neighbours. The numbers of nearest
! neighbours are stored in n(:). for example n(1)=4 means 4 nearest neighbours
!!! first neighb

       if (size(J_ij,2).ne.1) then
        lay=i_m
       else
        lay=1
       endif

      do i=1,param_N_Nneigh
       do j=1,indexNN(i,k)

        v_x=tableNN(1,avant+j,i_x,i_y,i_z,i_m)
        v_y=tableNN(2,avant+j,i_x,i_y,i_z,i_m)
        v_z=tableNN(3,avant+j,i_x,i_y,i_z,i_m)
        v_m=tableNN(4,avant+j,i_x,i_y,i_z,i_m)

       B_int=B_int+J_ij(i,lay)*dble(masque(avant+j+1,i_x,i_y,i_z)) &
        *Spin(X:Z,v_x,v_y,v_z,v_m)

       enddo
       avant=avant+indexNN(i,k)
      enddo

      if (size(indexNN,2).ne.1) then
      avant=sum(indexNN(:,1))
      k=2
      do i=1,param_N_Nneigh_il
       do j=1,indexNN(i,k)

       v_x=tableNN(1,avant+j,i_x,i_y,i_z,i_m)
       v_y=tableNN(2,avant+j,i_x,i_y,i_z,i_m)
       v_z=tableNN(3,avant+j,i_x,i_y,i_z,i_m)
       v_m=tableNN(4,avant+j,i_x,i_y,i_z,i_m)

       B_int=B_int+J_il(i)*dble(masque(avant+j+1,i_x,i_y,i_z)) &
        *Spin(X:Z,v_x,v_y,v_z,v_m)
      enddo
       avant=avant+indexNN(i,k)
      enddo
      endif

      if (dabs(J_z(1)).gt.1.0d-8) then
      avant=sum(indexNN(:,1))+sum(indexNN(:,2))
      k=3
      do i=1,count(dabs(J_z)>1.0d-8)
       do j=1,indexNN(i,k)

       v_x=tableNN(1,avant+j,i_x,i_y,i_z,i_m)
       v_y=tableNN(2,avant+j,i_x,i_y,i_z,i_m)
       v_z=tableNN(3,avant+j,i_x,i_y,i_z,i_m)
       v_m=tableNN(4,avant+j,i_x,i_y,i_z,i_m)

       B_int=B_int+J_z(i)*dble(masque(avant+j+1,i_x,i_y,i_z)) &
        *Spin(X:Z,v_x,v_y,v_z,v_m)
       enddo
       avant=avant+indexNN(i,k)
      enddo
      endif

      Bexch=-2.0d0*B_int*c_ji
      end function Bexch
!DM field
      function Bdm(i_x,i_y,i_z,i_m,spin,shape_spin,indexNN,shape_index,masque,shape_masque,tableNN,shape_tableNN,c_DM,DM)
      use m_vector, only: cross
      implicit none
! external variable
      integer , intent(in) :: i_x,i_y,i_z,i_m
      integer, intent(in) :: shape_index(2),shape_spin(5),shape_tableNN(6),shape_masque(4)
      integer, intent(in) :: tableNN(shape_tableNN(1),shape_tableNN(2),shape_tableNN(3),shape_tableNN(4),shape_tableNN(5),shape_tableNN(6))
      integer, intent(in) :: masque(shape_masque(1),shape_masque(2),shape_masque(3),shape_masque(4))
      integer, intent(in) :: indexNN(shape_index(1),shape_index(2))
      real(kind=8), intent(in) :: spin(shape_spin(1),shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5))
      real(kind=8), intent(in) :: c_DM,DM(:,:,:)
! value of the function
      real(kind=8), dimension(3) :: Bdm
!internal variable
      integer :: avant,i,j
      real(kind=8) :: B_int(3),DM_vector(6,3,1)
! coordinate of the spin
      integer :: X,Y,Z
! position of the neighbors along x,y,z and motif
      integer :: v_x,v_y,v_z,v_m

      avant=0
      X=shape_spin(1)-2
      Y=shape_spin(1)-1
      Z=shape_spin(1)-0

      B_int=0.0d0

      do i=1,count(abs(DM(:,i_m,1))>1.0d-8)
       do j=1,indexNN(i,1)

        v_x=tableNN(1,avant+j,i_x,i_y,i_z,i_m)
        v_y=tableNN(2,avant+j,i_x,i_y,i_z,i_m)
        v_z=tableNN(3,avant+j,i_x,i_y,i_z,i_m)
        v_m=tableNN(4,avant+j,i_x,i_y,i_z,i_m)

        B_int=B_int-DM(i,i_m,1)*cross(DM_vector(avant+j,:,i_m),Spin(X:Z,v_x,v_y,v_z,v_m)) &
         *dble(masque(avant+j+1,i_x,i_y,i_z))

       enddo
       avant=avant+indexNN(i,1)
      enddo

      Bdm=-2.0d0*B_int*c_DM

      end function Bdm
!!! Anisotropy field
      function Bani(i_x,i_y,i_z,i_m,spin,shape_spin,c_ani,D_ani)
      use m_vector, only : norm
      implicit none
! external variable
      integer , intent(in) :: i_x,i_y,i_z,i_m
      integer, intent(in) :: shape_spin(5)
      real(kind=8), intent(in) :: spin(shape_spin(1),shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5))
      real(kind=8) , intent(in) :: c_ani,D_ani(:,:,:)
! value of the function
      real(kind=8), dimension(3) :: Bani
! coordinate of the spin
      integer :: X,Y,Z
! internal variable
      real(kind=8) , dimension(3) :: B_int

      B_int=0.0d0
      X=shape_spin(1)-2
      Y=shape_spin(1)-1
      Z=shape_spin(1)-0

      B_int=B_int+2.0d0*(spin(X,i_x,i_y,i_z,i_m)*D_ani(1,1,1)+spin(Y,i_x,i_y,i_z,i_m)* &
      D_ani(2,1,1)+spin(Z,i_x,i_y,i_z,i_m)*D_ani(3,1,1))

      Bani=-c_ani*B_int
      end function Bani
! biquadratic field
      function Bbiqd(i_x,i_y,i_z,i_m,spin,shape_spin,indexNN,shape_index,masque,shape_masque,tableNN,shape_tableNN,c_JB,J_B)
      implicit none
! external variable
      integer , intent(in) :: i_x,i_y,i_z,i_m
      integer, intent(in) :: shape_index(2),shape_spin(5),shape_tableNN(6),shape_masque(4)
      integer, intent(in) :: tableNN(shape_tableNN(1),shape_tableNN(2),shape_tableNN(3),shape_tableNN(4),shape_tableNN(5),shape_tableNN(6))
      integer, intent(in) :: masque(shape_masque(1),shape_masque(2),shape_masque(3),shape_masque(4))
      integer, intent(in) :: indexNN(shape_index(1),shape_index(2))
      real(kind=8), intent(in) :: spin(shape_spin(1),shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5))
      real(kind=8), intent(in) :: c_JB,J_B
! value of the function
      real(kind=8), dimension(3) :: Bbiqd
! internal variable
      real(kind=8), dimension(3) :: B_int
! coordinate of the spin
      integer :: X,Y,Z
      integer :: j
! position of the neighbors along x,y,z and motif
      integer :: v_x,v_y,v_z,v_m

      B_int=0.0d0
      X=shape_spin(1)-2
      Y=shape_spin(1)-1
      Z=shape_spin(1)-0

! the choosen spin is the spin i_s. It has nn nearest neighbours. The numbers of nearest
! neighbours are stored in n(:). for example n(1)=4 means 4 nearest neighbours
!!! first neighb

      do j=1,indexNN(1,1)

        v_x=tableNN(1,j,i_x,i_y,i_z,i_m)
        v_y=tableNN(2,j,i_x,i_y,i_z,i_m)
        v_z=tableNN(3,j,i_x,i_y,i_z,i_m)
        v_m=tableNN(4,j,i_x,i_y,i_z,i_m)

       B_int=B_int+2.0d0*(Spin(X,i_x,i_y,i_z,i_m)*Spin(X,v_x,v_y,v_z,v_m)+ &
        Spin(Y,i_x,i_y,i_z,i_m)*Spin(Y,v_x,v_y,v_z,v_m)+ &
        Spin(Z,i_x,i_y,i_z,i_m)*Spin(Z,v_x,v_y,v_z,v_m)) &
        *Spin(X:Z,v_x,v_y,v_z,v_m)*dble(masque(j+1,i_x,i_y,i_z))
      enddo

      Bbiqd=-2.0d0*c_JB*J_B*B_int
      end function Bbiqd
!4 spin field
! the position of the corners are stored as a function of the lattice vectors i.e. the position
! ix,iy,iz. When corners is corners(k,1:3)=(/-1,1,0/), it meens ix-1, iy+1, iz+0
! ipuv is the corner along the diagonal u+v for k=3,6...
! ipu is along the first vector for k=1,4...
! ipv along the second k=2,5...
      function Bfour(i_x,i_y,i_z,i_m,spin,shape_spin,masque,shape_masque,my_lattice,c_Ki,K_1)
      use m_sym_utils, only : corners
      use m_derived_types
      implicit none
! external variable
      integer , intent(in) :: i_x,i_y,i_z,i_m
      integer, intent(in) :: shape_spin(5),shape_masque(4)
      integer, intent(in) :: masque(shape_masque(1),shape_masque(2),shape_masque(3),shape_masque(4))
      real(kind=8), intent(in) :: spin(shape_spin(1),shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5))
      real(kind=8), intent(in) :: c_Ki,K_1
      type(lattice), intent(in) :: my_lattice
! value of the function
      real(kind=8), dimension(3) :: Bfour
! coordinate of the spin
      integer :: X,Y,Z
! internal variable
      real(kind=8), dimension(3) :: B_int,B_local
      integer :: k,ipu(3),ipv(3),ipuv(3)
      logical :: Periodic_log(3)

      Periodic_log=my_lattice%boundary
      B_int=0.0d0
      X=shape_spin(1)-2
      Y=shape_spin(1)-1
      Z=shape_spin(1)-0

      k=1
      do while (k.lt.size(corners,1))

       B_local=0.0d0
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

       B_local=B_local+&
        Spin(X:Z,ipu(1),ipu(2),ipu(3),i_m)*& !ipu
         (Spin(X,ipv(1),ipv(2),ipv(3),i_m)*Spin(X,ipuv(1),ipuv(2),ipuv(3),i_m)+Spin(Y,ipv(1),ipv(2),ipv(3),i_m)* &
         Spin(Y,ipuv(1),ipuv(2),ipuv(3),i_m)+Spin(Z,ipv(1),ipv(2),ipv(3),i_m)*Spin(Z,ipuv(1),ipuv(2),ipuv(3),i_m))

       B_local=B_local+Spin(X:Z,ipv(1),ipv(2),ipv(3),i_m)*&
         (Spin(X,ipu(1),ipu(2),ipu(3),i_m)*Spin(X,ipuv(1),ipuv(2),ipuv(3),i_m)+Spin(Y,ipu(1),ipu(2),ipu(3),i_m)* &
         Spin(Y,ipuv(1),ipuv(2),ipuv(3),i_m)+Spin(Z,ipu(1),ipu(2),ipu(3),i_m)*Spin(Z,ipuv(1),ipuv(2),ipuv(3),i_m))

       B_local=B_local-Spin(X:Z,ipuv(1),ipuv(2),ipuv(3),i_m)*&
         (Spin(X,ipu(1),ipu(2),ipu(3),i_m)*Spin(X,ipv(1),ipv(2),ipv(3),i_m)+Spin(Y,ipu(1),ipu(2),ipu(3),i_m)* &
         Spin(Y,ipv(1),ipv(2),ipv(3),i_m)+Spin(Z,ipu(1),ipu(2),ipu(3),i_m)*Spin(Z,ipv(1),ipv(2),ipv(3),i_m))

       B_local=B_local*dble(masque(1,i_x,i_y,i_z)*masque(1,ipu(1),ipu(2),ipu(3))*&
       masque(1,ipv(1),ipv(2),ipv(3))*masque(1,ipuv(1),ipuv(2),ipuv(3)))
       B_int=B_int+B_local
       k=k+3
      enddo

      Bfour=-4.0d0*c_Ki*K_1*B_int
      end function Bfour
!Zeeman field
      function BZ(i_x,i_y,i_z,i_m,spin,shape_spin,h_ext)
      use m_constants, only : mu_B
      implicit none
! external variable
      integer , intent(in) :: i_x,i_y,i_z,i_m
      integer, intent(in) :: shape_spin(5)
      real(kind=8), intent(in) :: spin(shape_spin(1),shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5))
      real(kind=8), intent(in) :: h_ext(3)
! value of the function
      real(kind=8), dimension(3) :: BZ

      BZ=mu_B*H_ext*sqrt(Spin(1,i_x,i_y,i_z,i_m)**2+Spin(2,i_x,i_y,i_z,i_m)**2+Spin(3,i_x,i_y,i_z,i_m)**2)

      end function BZ
!!!! electric field
      function Efield(i_x,i_y,i_z,i_m,spin,shape_spin,indexNN,shape_index,masque,shape_masque,tableNN,shape_tableNN,c_ji,J_ij)
      use m_efield, only : me,efield_jij
      implicit none
! external variable
      integer , intent(in) :: i_x,i_y,i_z,i_m
      integer, intent(in) :: shape_index(2),shape_spin(5),shape_tableNN(6),shape_masque(4)
      integer, intent(in) :: tableNN(shape_tableNN(1),shape_tableNN(2),shape_tableNN(3),shape_tableNN(4),shape_tableNN(5),shape_tableNN(6))
      integer, intent(in) :: masque(shape_masque(1),shape_masque(2),shape_masque(3),shape_masque(4))
      integer, intent(in) :: indexNN(shape_index(1),shape_index(2))
      real(kind=8), intent(in) :: spin(shape_spin(1),shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5))
      real(kind=8), intent(in) :: c_ji,J_ij(:,:,:)
! value of the function
      real(kind=8), dimension(3) :: Efield
! internal variable
      integer :: i,j,avant
      real(kind=8), dimension(3) :: E_int
! coordinate of the spin
      integer :: X,Y,Z
! position of the neighbors along x,y,z and motif
      integer :: v_x,v_y,v_z,v_m

      E_int=0.0d0
      avant=0
      Efield=0.0d0
      X=shape_spin(1)-2
      Y=shape_spin(1)-1
      Z=shape_spin(1)-0

      if (sum(dabs(me)).lt.1.0d-8) return
! the choosen spin is the spin i_s. It has nn nearest neighbours. The numbers of nearest
! neighbours are stored in n(:). for example n(1)=4 means 4 nearest neighbours
!!! first neighb

      do i=1,count(dabs(J_ij)>1.0d-7)
       do j=1,indexNN(i,1)

        v_x=tableNN(1,avant+j,i_x,i_y,i_z,i_m)
        v_y=tableNN(2,avant+j,i_x,i_y,i_z,i_m)
        v_z=tableNN(3,avant+j,i_x,i_y,i_z,i_m)
        v_m=tableNN(4,avant+j,i_x,i_y,i_z,i_m)

       E_int=E_int+me(i)*Efield_Jij(i_x,i_y,i_z)*Spin(X:Z,v_x,v_y,v_z,v_m)* &
        dble(masque(avant+j+1,i_x,i_y,i_z))
       enddo
       avant=avant+indexNN(i,1)
      enddo

      Efield=-2.0d0*E_int*c_ji
      end function Efield

#ifdef CPP_BRUTDIP
! Dipole Dipole interaction
      function Bdip(i_x,i_y,i_z,i_m,spin,shape_spin)
      use m_constants, only : pi
      use m_vector, only : norm
      implicit none
! value of the function
      real(kind=8) :: Bdip(3)
! external variable
      integer , intent(in) :: i_x,i_y,i_z,i_m
      integer, intent(in) :: shape_spin(5)
      real(kind=8), intent(in) :: spin(shape_spin(1),shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5))
! internal variable
      integer :: j_x,j_y,j_z,j_m,nmag
      real(kind=8) :: rc(3),step(3),ss
      real(kind=8), parameter :: alpha=6.74582d-7

      Bdip=0.0d0
      nmag=shape_spin(5)

! the choosen spin is the spin i_s. It has nn nearest neighbours. The numbers of nearest
! neighbours are stored in n(:). for example n(1)=4 means 4 nearest neighbours
!!! first neighb
       do j_m=1,shape_spin(5)
        do j_z=1,shape_spin(4)
         do j_y=1,shape_spin(3)
          do j_x=1,shape_spin(2)

          rc=spin(1:3,j_x,j_y,j_z,j_m)-spin(1:3,i_x,i_y,i_z,i_m)

        ss=norm(rc)
        if (ss.lt.1.0d-3) cycle
        rc=rc/ss

        step=rc*dot_product(rc,spin(4:6,j_x,j_y,j_z,j_m))

        Bdip=Bdip+(spin(4:6,j_x,j_y,j_z,j_m)-3.0d0*step)/ss**3

          enddo
         enddo
        enddo
       enddo

      Bdip=-2.0d0*Bdip/pi(4.0d0)*0.5d0*alpha

      end function Bdip
#else

! Dipole Dipole interaction with FFT
      function Bdip(i_x,i_y,i_z,i_m,spin,shape_spin)
      use m_setup_dipole, only : mmatrix,mcomplex,hreal,hcomplex,Nfftx,Nffty,Nfftz,ntensor &
     & ,rtrans,ctrans,planrtoc,planctor
      use m_fft
      use m_constants, only : mu_B
      implicit none
! value of the function
      real(kind=8) :: Bdip(3)
! external variable
      integer , intent(in) :: i_x,i_y,i_z,i_m
      integer, intent(in) :: shape_spin(5),shape_masque(4)
      integer, intent(in) :: masque(shape_masque(1),shape_masque(2),shape_masque(3),shape_masque(4))
      real(kind=8), intent(in) :: spin(shape_spin(1),shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5))
! internal variable
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

      Bdip=hreal(:,i_x,i_y,i_z)*mu_B

      end function Bdip
#endif
      end module m_fieldeff
