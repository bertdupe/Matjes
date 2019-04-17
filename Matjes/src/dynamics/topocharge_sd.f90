module m_topo_sd
      use m_derived_types
      interface topocharge_sd
       module procedure topo_2D
      end interface topocharge_sd
      contains

      subroutine topo_2D(j,spin,my_lattice)
      use m_vector, only : cross,norm,area
      use m_constants, only : pi
      implicit none
      type(lattice), intent(in) :: my_lattice
      integer, intent(in) :: j
      real(kind=8), intent(in) :: spin(:,:,:,:)

!internal
      real(kind=8), allocatable :: map_toto(:,:)
      real(kind=8), allocatable :: map_vort(:,:,:)
      integer :: i_x,i_y,ipu,ipv,k,i_m,ipm
      integer :: i,dim_lat(3),shape_spin(4)
      character(len=30) :: fname,toto
      real(kind=8) :: surface,net(3,3)
      logical :: Periodic_log(3)


      dim_lat=my_lattice%dim_lat
      net=my_lattice%areal
      allocate(map_toto(dim_lat(1),dim_lat(2)),map_vort(3,dim_lat(1),dim_lat(2)))
      surface=norm(cross(net(1,:),net(2,:)))
      map_toto=0.0d0
      map_vort=0.0d0
      shape_spin=shape(spin)

      if (shape_spin(4).eq.1) then
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

       map_toto(i_x,i_y)=area(spin(:,i_x,i_y,1),spin(:,ipu,i_y,1),spin(:,ipu,ipv,1))+ &
          area(spin(:,i_x,i_y,1),spin(:,ipu,ipv,1),spin(:,i_x,ipv,1))

! Calcul de la vorticite

       map_vort(:,i_x,i_y)=cross(Spin(:,i_x,i_y,1),Spin(:,ipu,i_y,1))+ &
        cross(Spin(:,ipu,i_y,1),Spin(:,ipu,ipv,1))+ &
        cross(Spin(:,ipu,ipv,1),Spin(:,i_x,ipv,1))+ &
        cross(Spin(:,i_x,ipv,1),Spin(:,i_x,i_y,1))

        enddo
       enddo
#ifdef CPP_OPENMP
!$OMP end parallel do
#endif
      else

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
        map_toto(i_x,i_y)=map_toto(i_x,i_y)+area(spin(:,i_x,i_y,i_m),spin(:,ipu,i_y,i_m),spin(:,ipu,ipv,i_m))+ &
         area(spin(:,i_x,i_y,i_m),spin(:,ipu,ipv,i_m),spin(:,i_x,ipv,i_m))

!!!!!!!!!!!!!!!!!!
! bottom surface
!!!!!!!!!!!!!!!!!!
        i_m=1

        map_toto(i_x,i_y)=map_toto(i_x,i_y)+area(spin(:,i_x,i_y,i_m),spin(:,ipu,i_y,i_m),spin(:,ipu,ipv,i_m))+ &
         area(spin(:,i_x,i_y,i_m),spin(:,ipu,ipv,i_m),spin(:,i_x,ipv,i_m))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! unitl the end, i_m=1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        i_m=1
!!!!!!!!!!!!!!!!!!
! (0-10) surface (make a drawing...) the rotation is from u to w
!!!!!!!!!!!!!!!!!!
        ipm=2

        map_toto(i_x,i_y)=map_toto(i_x,i_y)+area(spin(:,i_x,i_y,i_m),spin(:,ipu,i_y,i_m),spin(:,ipu,i_y,ipm))+ &
         area(spin(:,i_x,i_y,i_m),spin(:,ipu,i_y,ipm),spin(:,i_x,i_y,ipm))
!!!!!!!!!!!!!!!!!!
! (100) surface (make a drawing...) the rotation is from v to w
!!!!!!!!!!!!!!!!!!
        map_toto(i_x,i_y)=map_toto(i_x,i_y)+area(spin(:,ipu,i_y,i_m),spin(:,ipu,ipv,i_m),spin(:,ipu,ipv,ipm))+ &
         area(spin(:,ipu,i_y,i_m),spin(:,ipu,ipv,ipm),spin(:,ipu,i_y,ipm))
!!!!!!!!!!!!!!!!!!
! (010) surface (make a drawing...) the rotation is from -u to w
!!!!!!!!!!!!!!!!!!
        map_toto(i_x,i_y)=map_toto(i_x,i_y)+area(spin(:,ipu,ipv,i_m),spin(:,i_x,ipv,i_m),spin(:,i_x,ipv,ipm))+ &
         area(spin(:,ipu,ipv,i_m),spin(:,i_x,ipv,ipm),spin(:,ipu,ipv,ipm))
!!!!!!!!!!!!!!!!!!
! (-100) surface (make a drawing...) the rotation is from w to v
!!!!!!!!!!!!!!!!!!
        map_toto(i_x,i_y)=map_toto(i_x,i_y)+area(spin(:,i_x,i_y,i_m),spin(:,i_x,i_y,ipm),spin(:,i_x,ipv,ipm))+ &
         area(spin(:,i_x,i_y,i_m),spin(:,i_x,ipv,ipm),spin(:,i_x,ipv,i_m))

! Calcul de la vorticite

       map_vort(:,i_x,i_y)=cross(Spin(4:6,i_x,i_y,1),Spin(4:6,ipu,i_y,1))+ &
        cross(Spin(4:6,ipu,i_y,1),Spin(4:6,ipu,ipv,1))+ &
        cross(Spin(4:6,ipu,ipv,1),Spin(4:6,i_x,ipv,1))+ &
        cross(Spin(4:6,i_x,ipv,1),Spin(4:6,i_x,i_y,1))

        enddo
       enddo
#ifdef CPP_OPENMP
!$OMP end parallel do
#endif

      endif

      write(fname,'(I12)') j
      toto=trim(adjustl(fname))
      write(fname,'(a,18a,a)')'topomap',(toto(i:i),i=1, &
       len_trim(toto)),'.dat'
      OPEN(70,FILE=fname,action='write',status='unknown',form='formatted')

      do i_y=1,dim_lat(2)
       do i_x=1,dim_lat(1)
        Write(70,'(6(f14.8,2x))') (Spin(k,i_x,i_y,1),k=1,2), map_vort(:,i_x,i_y)/3.0d0/dsqrt(3.0d0), &
     &   map_toto(i_x,i_y)/pi(4.0d0)
       enddo
      enddo
      close(70)

      end subroutine topo_2D

end module m_topo_sd
