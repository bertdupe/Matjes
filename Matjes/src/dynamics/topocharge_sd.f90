      module m_topo_sd
      interface topocharge_sd
       module procedure topo_2D
       module procedure topo_2D_SL
       module procedure topo_1D
       module procedure topo_3D
      end interface topocharge_sd
      contains

      subroutine topo_2D(j,spin)
      use m_rw_lattice, only : dim_lat,net
      use m_vector, only : cross,norm,area
      use m_constants, only : pi
      use m_parameters, only : Periodic_log
      implicit none
      integer, intent(in) :: j
      real(kind=8), intent(in) :: spin(:,:,:)
      real(kind=8), dimension(dim_lat(1),dim_lat(2)) :: map_toto
      real(kind=8), dimension(3,dim_lat(1),dim_lat(2)) :: map_vort
      integer :: i_x,i_y,ipu,ipv,k
      integer :: iomp,i
      character(len=30) :: fname,toto
      real(kind=8) :: surface

      surface=norm(cross(net(1,:),net(2,:)))
      map_toto=0.0d0
      map_vort=0.0d0

#ifdef CPP_OPENMP
!$OMP parallel do private(i_x,i_y,ipu,ipv) default(shared)
#endif
       do i_x=1,dim_lat(1)
        do i_y=1,dim_lat(2)

        if (all(Periodic_log)) then
         ipu=mod(i_x-1+1+dim_lat(1),dim_lat(1))+1
         ipv=mod(i_y+1+dim_lat(2)-1,dim_lat(2))+1
        else
         if (i_x.eq.dim_lat(1)) then
          ipu=i_x
          else
          ipu=i_x+1
         endif
         if (i_y.eq.dim_lat(2)) then
          ipv=i_y
          else
          ipv=i_y+1
         endif
        endif

! charge euler

       map_toto(i_x,i_y)=area(spin(4:6,i_x,i_y),spin(4:6,ipu,i_y),spin(4:6,ipu,ipv))+ &
          area(spin(4:6,i_x,i_y),spin(4:6,ipu,ipv),spin(4:6,i_x,ipv))

! Calcul de la vorticite

       map_vort(:,i_x,i_y)=cross(Spin(4:6,i_x,i_y),Spin(4:6,ipu,i_y))+ &
        cross(Spin(4:6,ipu,i_y),Spin(4:6,ipu,ipv))+ &
        cross(Spin(4:6,ipu,ipv),Spin(4:6,i_x,ipv))+ &
        cross(Spin(4:6,i_x,ipv),Spin(4:6,i_x,i_y))

        enddo
       enddo
#ifdef CPP_OPENMP
!$OMP end parallel do
#endif
      write(fname,'(I12)') j
      toto=trim(adjustl(fname))
      write(fname,'(a,18a,a)')'topomap',(toto(i:i),i=1, &
       len_trim(toto)),'.dat'
      OPEN(70,FILE=fname,action='write',status='unknown',form='formatted')

      do i_x=1,dim_lat(1)
       do i_y=1,dim_lat(2)
        Write(70,'(6(f14.8,2x))') (Spin(k,i_x,i_y),k=1,2), map_vort(:,i_x,i_y)/3.0d0/dsqrt(3.0d0), &
     &   map_toto(i_x,i_y)/pi(4.0d0)
       enddo
      enddo
      close(70)

      end subroutine topo_2D

      subroutine topo_2D_SL(j,spins)
      use m_rw_lattice, only : dim_lat
      use m_vector, only : cross,norm,area
      use m_constants, only : pi
      use m_parameters, only : Periodic_log
      implicit none
      integer, intent(in) :: j
      real(kind=8), intent(in) :: spins(:,:,:,:)
      real(kind=8), dimension(dim_lat(1),dim_lat(2)) :: map_toto
      real(kind=8), dimension(3,dim_lat(1),dim_lat(2)) :: map_vort
      integer :: i_x,i_y,ipu,ipv,k,i_m,ipm
      integer :: iomp,i
      character(len=30) :: fname,toto

      map_toto=0.0d0
      map_vort=0.0d0

#ifdef CPP_OPENMP
!$OMP parallel do private(i_x,i_y,i_m,ipm,ipu,ipv) default(shared)
#endif
       do i_x=1,dim_lat(1)
        do i_y=1,dim_lat(2)

        if (all(Periodic_log)) then
         ipu=mod(i_x-1+1+dim_lat(1),dim_lat(1))+1
         ipv=mod(i_y+1+dim_lat(2)-1,dim_lat(2))+1
        else
         if (i_x.eq.dim_lat(1)) then
          ipu=i_x
          else
          ipu=i_x+1
         endif
         if (i_y.eq.dim_lat(2)) then
          ipv=i_y
          else
          ipv=i_y+1
         endif
        endif

! charge euler
        i_m=2
!!!!!!!!!!!!!!!!!!
! above surface
!!!!!!!!!!!!!!!!!!
        map_toto(i_x,i_y)=map_toto(i_x,i_y)+area(spins(:,i_x,i_y,i_m),spins(:,ipu,i_y,i_m),spins(:,ipu,ipv,i_m))+ &
         area(spins(:,i_x,i_y,i_m),spins(:,ipu,ipv,i_m),spins(:,i_x,ipv,i_m))

!!!!!!!!!!!!!!!!!!
! bottom surface
!!!!!!!!!!!!!!!!!!
        i_m=1

        map_toto(i_x,i_y)=map_toto(i_x,i_y)+area(spins(:,i_x,i_y,i_m),spins(:,ipu,i_y,i_m),spins(:,ipu,ipv,i_m))+ &
         area(spins(:,i_x,i_y,i_m),spins(:,ipu,ipv,i_m),spins(:,i_x,ipv,i_m))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! unitl the end, i_m=1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        i_m=1
!!!!!!!!!!!!!!!!!!
! (0-10) surface (make a drawing...) the rotation is from u to w
!!!!!!!!!!!!!!!!!!
        ipm=2

        map_toto(i_x,i_y)=map_toto(i_x,i_y)+area(spins(:,i_x,i_y,i_m),spins(:,ipu,i_y,i_m),spins(:,ipu,i_y,ipm))+ &
         area(spins(:,i_x,i_y,i_m),spins(:,ipu,i_y,ipm),spins(:,i_x,i_y,ipm))
!!!!!!!!!!!!!!!!!!
! (100) surface (make a drawing...) the rotation is from v to w
!!!!!!!!!!!!!!!!!!
        map_toto(i_x,i_y)=map_toto(i_x,i_y)+area(spins(:,ipu,i_y,i_m),spins(:,ipu,ipv,i_m),spins(:,ipu,ipv,ipm))+ &
         area(spins(:,ipu,i_y,i_m),spins(:,ipu,ipv,ipm),spins(:,ipu,i_y,ipm))
!!!!!!!!!!!!!!!!!!
! (010) surface (make a drawing...) the rotation is from -u to w
!!!!!!!!!!!!!!!!!!
        map_toto(i_x,i_y)=map_toto(i_x,i_y)+area(spins(:,ipu,ipv,i_m),spins(:,i_x,ipv,i_m),spins(:,i_x,ipv,ipm))+ &
         area(spins(:,ipu,ipv,i_m),spins(:,i_x,ipv,ipm),spins(:,ipu,ipv,ipm))
!!!!!!!!!!!!!!!!!!
! (-100) surface (make a drawing...) the rotation is from w to v
!!!!!!!!!!!!!!!!!!
        map_toto(i_x,i_y)=map_toto(i_x,i_y)+area(spins(:,i_x,i_y,i_m),spins(:,i_x,i_y,ipm),spins(:,i_x,ipv,ipm))+ &
         area(spins(:,i_x,i_y,i_m),spins(:,i_x,ipv,ipm),spins(:,i_x,ipv,i_m))

! Calcul de la vorticite

       map_vort(:,i_x,i_y)=cross(Spins(4:6,i_x,i_y,1),Spins(4:6,ipu,i_y,1))+ &
        cross(Spins(4:6,ipu,i_y,1),Spins(4:6,ipu,ipv,1))+ &
        cross(Spins(4:6,ipu,ipv,1),Spins(4:6,i_x,ipv,1))+ &
        cross(Spins(4:6,i_x,ipv,1),Spins(4:6,i_x,i_y,1))

        enddo
       enddo
#ifdef CPP_OPENMP
!$OMP end parallel do
#endif

      write(fname,'(I6)') j
      toto=trim(adjustl(fname))
      write(fname,'(a,18a,a)')'topomap',(toto(i:i),i=1, &
       len_trim(toto)),'.dat'
      OPEN(70,FILE=fname,action='write',status='unknown', &
       position='append',form='formatted')

      do i_x=1,dim_lat(1)
       do i_y=1,dim_lat(2)
        Write(70,'(6(f14.8,2x))') (Spins(k,i_x,i_y,1),k=1,2), map_vort(:,i_x,i_y)/3.0d0/dsqrt(3.0d0), &
         map_toto(i_x,i_y)/pi(4.0d0)
       enddo
      enddo
      close(70)

      end subroutine topo_2D_SL
!-----------------------------
! todo
      subroutine topo_1D(j,spins)
      use m_rw_lattice, only : dim_lat
      use m_vector, only : cross, norm
      use m_constants, only : pi
      implicit none
      integer, intent(in) :: j
      real(kind=8), intent(in) :: spins(:,:)
      real(kind=8), dimension(dim_lat(1)) :: map_toto
      real(kind=8), dimension(3,dim_lat(1)) :: map_vort
      integer :: i_x,i_y,i_z,ipu,ipv,k,i_m,ipm
      integer :: iomp,i
      character(len=30) :: fname,toto

      map_toto=0.0d0
      map_vort=0.0d0

      end subroutine topo_1D

      subroutine topo_3D(j,spins)
      use m_rw_lattice, only : dim_lat
      use m_vector, only : cross, norm
      use m_constants, only : pi
      implicit none
      integer, intent(in) :: j
      real(kind=8), intent(in) :: spins(:,:,:,:,:)
      real(kind=8), dimension(dim_lat(1),dim_lat(2),dim_lat(3)) :: map_toto
      real(kind=8), dimension(3,dim_lat(1),dim_lat(2),dim_lat(3)) :: map_vort
      integer :: i_x,i_y,i_z,ipu,ipv,k,i_m,ipm
      integer :: iomp,i
      character(len=30) :: fname,toto

      map_toto=0.0d0
      map_vort=0.0d0

      end subroutine topo_3D

      end module m_topo_sd
