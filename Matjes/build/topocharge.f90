      module m_topocharge
      interface topo_map
       module procedure topo_2D
       module procedure topo_2D_SL
       module procedure topo_1D
      end interface topo_map
      contains

      subroutine topo_2D(spins,shape_spin,signature,i_plot,Periodic_log)
      use m_vector, only : cross, area
      use m_topoplot
#ifdef CPP_MPI
      use m_parameters, only : i_ghost,n_ghost
      use m_mpi
      use m_mpi_prop, only : MPI_COMM_BOX
      use m_make_box, only : Xstart,Xstop,Ystart,Ystop
#endif
      implicit none
      real(kind=8), intent(in) :: spins(:,:,:),signature
      integer, intent(in) :: shape_spin(:)
      logical, intent(in) :: i_plot,Periodic_log(:)
!internal
      integer :: i,j,ipu,ipv
      real(kind=8) :: map_eul_int(shape_spin(2),shape_spin(3)),map_vort_int(3,shape_spin(2),shape_spin(3))
      real(kind=8) :: map_eul(shape_spin(2),shape_spin(3)),map_vort(3,shape_spin(2),shape_spin(3))
#ifndef CPP_MPI
      integer :: Xstart,Xstop,Ystart,Ystop

      Xstart=1
      Xstop=shape_spin(2)
      Ystart=1
      Ystop=shape_spin(3)
#endif

      map_eul_int=0.0d0
      map_vort_int=0.0d0

#ifdef CPP_OPENMP
!$OMP parallel private(i,j) default(shared)
#endif


       do j=Ystart,Ystop
        do i=Xstart,Xstop

        if (Periodic_log(1).and.Periodic_log(2)) then
         ipu=mod(i+shape_spin(2),shape_spin(2))+1
         ipv=mod(j+shape_spin(3),shape_spin(3))+1
        else
         if (i.eq.shape_spin(2)) then
          ipu=i
          else
          ipu=i+1
         endif
         if (j.eq.shape_spin(3)) then
          ipv=j
          else
          ipv=j+1
         endif
        endif

       map_eul_int(i,j)=area(Spins(:,i,j),Spins(:,ipu,j),Spins(:,ipu,ipv))+ &
        area(Spins(:,i,j),Spins(:,ipu,ipv),Spins(:,i,ipv))

       map_vort_int(:,i,j)=cross(Spins(:,i,j),Spins(:,ipu,j))+ &
       cross(Spins(:,ipu,j),Spins(:,ipu,ipv))+cross(Spins(:,ipu,ipv),Spins(:,i,ipv))+ &
       cross(Spins(:,i,ipv),Spins(:,i,j))

        enddo
       enddo

#ifdef CPP_OPENMP
!$OMP end parallel
#endif

#ifdef CPP_MPI
      if (i_ghost) map_eul=allreduce(map_eul_int,(/Xstop-Xstart+1,Ystop-Ystart+1/),MPI_COMM_BOX)
      if (i_ghost) map_vort=allreduce(map_vort_int,(/Xstop-Xstart+1,Ystop-Ystart+1/),3,MPI_COMM_BOX)
#else
      map_eul=map_eul_int
      map_vort=map_vort_int
#endif

      if (i_plot) call topoplot(signature,map_eul,map_vort)

      end subroutine topo_2D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine topo_2D_SL(spins,shape_spin,signature,i_plot,Periodic_log)
      use m_vector, only : cross, area
      use m_topoplot
#ifdef CPP_MPI
      use m_parameters, only : i_ghost,n_ghost
      use m_mpi
      use m_mpi_prop, only : MPI_COMM_BOX
      use m_make_box, only : Xstart,Xstop,Ystart,Ystop
#endif
      implicit none
      real(kind=8), intent(in) :: spins(:,:,:,:),signature
      integer, intent(in) :: shape_spin(:)
      logical, intent(in) :: i_plot,Periodic_log(:)
!internal
      integer :: i_x,i_y,ipu,ipv,i_m,ipm
      real(kind=8) :: map_eul_int(shape_spin(2),shape_spin(3)),map_vort_int(3,shape_spin(2),shape_spin(3))
      real(kind=8) :: map_eul(shape_spin(2),shape_spin(3)),map_vort(3,shape_spin(2),shape_spin(3))
#ifndef CPP_MPI
      integer :: Xstart,Xstop,Ystart,Ystop

      Xstart=1
      Xstop=shape_spin(2)
      Ystart=1
      Ystop=shape_spin(3)
#endif

      map_eul_int=0.0d0
      map_vort_int=0.0d0


#ifdef CPP_OPENMP
!$OMP parallel private(i_x,i_y,i_m,ipm) default(shared)
#endif

      do i_y=Ystart,Ystop
       do i_x=Xstart,Xstop

        if (Periodic_log(1).and.Periodic_log(2)) then
         ipu=mod(i_x-1+1+shape_spin(2),shape_spin(2))+1
         ipv=mod(i_y+1+shape_spin(3)-1,shape_spin(3))+1
        else
         if (i_x.eq.shape_spin(2)) then
          ipu=i_x
          else
          ipu=i_x+1
         endif
         if (i_y.eq.shape_spin(3)) then
          ipv=i_y
          else
          ipv=i_y+1
         endif
        endif

        i_m=2
!!!!!!!!!!!!!!!!!!
! above surface
!!!!!!!!!!!!!!!!!!
        map_eul_int(i_x,i_y)=map_eul_int(i_x,i_y)+area(spins(:,i_x,i_y,i_m),spins(:,ipu,i_y,i_m),spins(:,ipu,ipv,i_m))+ &
         area(spins(:,i_x,i_y,i_m),spins(:,ipu,ipv,i_m),spins(:,i_x,ipv,i_m))

!!!!!!!!!!!!!!!!!!
! bottom surface
!!!!!!!!!!!!!!!!!!
        i_m=1

        map_eul_int(i_x,i_y)=map_eul_int(i_x,i_y)+area(spins(:,i_x,i_y,i_m),spins(:,ipu,i_y,i_m),spins(:,ipu,ipv,i_m))+ &
         area(spins(:,i_x,i_y,i_m),spins(:,ipu,ipv,i_m),spins(:,i_x,ipv,i_m))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! unitl the end, i_m=1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        i_m=1
!!!!!!!!!!!!!!!!!!
! (0-10) surface (make a drawing...) the rotation is from u to w
!!!!!!!!!!!!!!!!!!
        ipm=2

        map_eul_int(i_x,i_y)=map_eul_int(i_x,i_y)+area(spins(:,i_x,i_y,i_m),spins(:,ipu,i_y,i_m),spins(:,ipu,i_y,ipm))+ &
         area(spins(:,i_x,i_y,i_m),spins(:,ipu,i_y,ipm),spins(:,i_x,i_y,ipm))
!!!!!!!!!!!!!!!!!!
! (100) surface (make a drawing...) the rotation is from v to w
!!!!!!!!!!!!!!!!!!
        map_eul_int(i_x,i_y)=map_eul_int(i_x,i_y)+area(spins(:,ipu,i_y,i_m),spins(:,ipu,ipv,i_m),spins(:,ipu,ipv,ipm))+ &
         area(spins(:,ipu,i_y,i_m),spins(:,ipu,ipv,ipm),spins(:,ipu,i_y,ipm))
!!!!!!!!!!!!!!!!!!
! (010) surface (make a drawing...) the rotation is from -u to w
!!!!!!!!!!!!!!!!!!
        map_eul_int(i_x,i_y)=map_eul_int(i_x,i_y)+area(spins(:,ipu,ipv,i_m),spins(:,i_x,ipv,i_m),spins(:,i_x,ipv,ipm))+ &
         area(spins(:,ipu,ipv,i_m),spins(:,i_x,ipv,ipm),spins(:,ipu,ipv,ipm))
!!!!!!!!!!!!!!!!!!
! (-100) surface (make a drawing...) the rotation is from w to v
!!!!!!!!!!!!!!!!!!
        map_eul_int(i_x,i_y)=map_eul_int(i_x,i_y)+area(spins(:,i_x,i_y,i_m),spins(:,i_x,i_y,ipm),spins(:,i_x,ipv,ipm))+ &
         area(spins(:,i_x,i_y,i_m),spins(:,i_x,ipv,ipm),spins(:,i_x,ipv,i_m))

       enddo
      enddo
#ifdef CPP_OPENMP
!$OMP end parallel

!$OMP parallel private(i_x,i_y,i_m,ipm) default(shared)
#endif
      do i_m=1,shape_spin(5)
       do i_y=Ystart,Ystop
        do i_x=Xstart,Xstop

        if (Periodic_log(1).and.Periodic_log(2)) then
         ipu=mod(i_x-1+1+shape_spin(2),shape_spin(2))+1
         ipv=mod(i_y+1+shape_spin(3)-1,shape_spin(3))+1
        else
         if (i_x.eq.shape_spin(2)) then
          ipu=i_x
          else
          ipu=i_x+1
         endif
         if (i_y.eq.shape_spin(3)) then
          ipv=i_y
          else
          ipv=i_y+1
         endif
        endif

        map_vort_int(:,i_x,i_y)=map_vort_int(:,i_x,i_y)+cross(Spins(:,i_x,i_y,i_m),Spins(:,ipu,i_y,i_m))+ &
         cross(Spins(:,ipu,i_y,i_m),Spins(:,ipu,ipv,i_m))+cross(Spins(:,ipu,ipv,i_m),Spins(:,i_x,ipv,i_m))+ &
         cross(Spins(:,i_x,ipv,i_m),Spins(:,i_x,i_y,i_m))

        enddo
       enddo
      enddo

#ifdef CPP_OPENMP
!$OMP end parallel
#endif

#ifdef CPP_MPI
      if (i_ghost) map_eul=allreduce(map_eul_int,(/Xstop-Xstart+1,Ystop-Ystart+1/),MPI_COMM_BOX)
      if (i_ghost) map_vort=allreduce(map_vort_int,(/Xstop-Xstart+1,Ystop-Ystart+1/),3,MPI_COMM_BOX)
#else
      map_eul=map_eul_int
      map_vort=map_vort_int
#endif

      if (i_plot) call topoplot(signature,map_eul,map_vort)

      end subroutine topo_2D_SL

!!!!!!!!!!!!!!!!!!!!!!!!! 1D case
      subroutine topo_1D(spins,shape_spin,signature,i_plot,Periodic_log)
      use m_topoplot
      use m_vector, only : cross,norm
#ifdef CPP_MPI
      use m_parameters, only : i_ghost,n_ghost
      use m_mpi
      use m_mpi_prop, only : MPI_COMM_BOX
      use m_make_box, only : Xstart,Xstop,Ystart,Ystop
#endif
      implicit none
      real(kind=8), intent(in) :: spins(:,:),signature
      integer, intent(in) :: shape_spin(:)
      logical, intent(in) :: i_plot,Periodic_log(:)
!internal
      integer :: i,ipu
      real(kind=8) :: map_eul_int(shape_spin(2)),map_vort_int(3,shape_spin(2))
      real(kind=8) :: map_eul(shape_spin(2)),map_vort(3,shape_spin(2))
      real(kind=8) :: dumy
#ifndef CPP_MPI
      integer :: Xstart,Xstop

      Xstart=1
      Xstop=shape_spin(2)
#endif

      map_eul_int=0.0d0
      map_vort_int=0.0d0

#ifdef CPP_OPENMP
!$OMP parallel do default(shared)
#endif


       do i=Xstart,Xstop

        if (Periodic_log(1)) then
         ipu=mod(i+shape_spin(2),shape_spin(2))+1
        else
         if (i.eq.shape_spin(2)) then
          ipu=i
          else
          ipu=i+1
         endif
        endif

        dumy=acos(Spins(1,i)*Spins(1,ipu)+Spins(2,i)*Spins(2,ipu)+Spins(3,i)*Spins(3,ipu))
        map_eul_int(i)=dumy*(Spins(3,i)+Spins(3,ipu))*norm(Spins(:,i)+Spins(:,ipu))

        map_vort_int(:,i)=cross(Spins(:,i),Spins(:,ipu))

       enddo

#ifdef CPP_OPENMP
!$OMP end parallel
#endif

#ifdef CPP_MPI
      if (i_ghost) map_eul=allreduce(map_eul_int,Xstop-Xstart+1,MPI_COMM_BOX)
      if (i_ghost) map_vort=allreduce(map_vort_int,(/3,Xstop-Xstart+1/),MPI_COMM_BOX)
#else
      map_eul=map_eul_int
      map_vort=map_vort_int
#endif

      if (i_plot) call topoplot(signature,map_eul,map_vort)

      end subroutine topo_1D

      end module m_topocharge
