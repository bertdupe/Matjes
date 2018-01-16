      module m_qorien
      interface qorien
       module procedure qorien_2D
       module procedure qorien_2D_SL
      end interface qorien
      contains

      subroutine qorien_2D(spins,shape_spin,my_lattice)
      use m_vector, only : sorien
      use m_derived_types
#ifdef CPP_MPI
      use m_parameters, only : i_ghost,Periodic_log,n_ghost
      use m_mpi
      use m_make_box, only : Xstart,Xstop,Ystart,Ystop
      use m_mpi_prop, only : MPI_COMM_BOX
#endif
      implicit none
      real(kind=8), intent(in) :: spins(:,:,:)
      integer, intent(in) :: shape_spin(:)
      type(lattice), intent(in) :: my_lattice
!internal
      integer :: i,j,ipu,ipv
      real(kind=8) :: map_int(3,shape_spin(2),shape_spin(3))
      real(kind=8) :: map_qorien(3,shape_spin(2),shape_spin(3))
! coordinate of the spin
      integer :: X,Y,Z
      logical :: Periodic_log(3)
#ifndef CPP_MPI
      integer :: Xstart,Xstop,Ystart,Ystop

      Xstart=1
      Xstop=shape_spin(2)
      Ystart=1
      Ystop=shape_spin(3)
#endif

      map_int=0.0d0
      X=shape_spin(1)-3
      Y=shape_spin(1)-2
      Z=shape_spin(1)-1
      Periodic_log=my_lattice%boundary

#ifdef CPP_OPENMP
!$OMP parallel DO private(i,j) default(shared)
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

       map_int(:,i,j)=(sorien(Spins(X:Z,i,j),Spins(X:Z,ipu,j),Spins(X:Z,ipu,ipv))+ &
        sorien(Spins(X:Z,i,j),Spins(X:Z,ipu,ipv),Spins(X:Z,i,ipv)))/sqrt(2.0d0)

        enddo
       enddo

#ifdef CPP_OPENMP
!$OMP end parallel do
#endif

#ifdef CPP_MPI
      if (i_ghost) then
       map_qorien=allreduce(map_int,shape_spin(2:3),3,MPI_COMM_BOX)
       else
       map_qorien=map_int
      endif
#else
       map_qorien=map_int
#endif

      end subroutine qorien_2D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine qorien_2D_SL(spins,shape_spin,my_lattice)
      use m_vector, only : sorien
      use m_derived_types
#ifdef CPP_MPI
      use m_parameters, only : i_ghost,Periodic_log,n_ghost
      use m_mpi
      use m_make_box, only : Xstart,Xstop,Ystart,Ystop
      use m_mpi_prop, only : MPI_COMM_BOX
#endif
      implicit none
      real(kind=8), intent(in) :: spins(:,:,:,:)
      integer, intent(in) :: shape_spin(:)
      type(lattice), intent(in) :: my_lattice
!internal
      integer :: i_x,i_y,ipu,ipv,i_m,ipm
      real(kind=8) :: map_int(3,shape_spin(2),shape_spin(3))
      real(kind=8) :: map_qorien(3,shape_spin(2),shape_spin(3))
! coordinate of the spin
      integer :: X,Y,Z
#ifndef CPP_MPI
      integer :: Xstart,Xstop,Ystart,Ystop
      logical :: Periodic_log(3)

      Periodic_log=my_lattice%boundary
      Xstart=1
      Xstop=shape_spin(2)
      Ystart=1
      Ystop=shape_spin(3)
#endif

      map_int=0.0d0
      X=shape_spin(1)-3
      Y=shape_spin(1)-2
      Z=shape_spin(1)-1
#ifdef CPP_OPENMP
!$OMP parallel DO private(i_x,i_y) default(shared)
#endif

      do i_y=Ystart,Ystop
       do i_x=Xstop,Xstart

        if (Periodic_log(1).and.Periodic_log(2)) then
         ipu=mod(i_x+shape_spin(2),shape_spin(2))+1
         ipv=mod(i_y+shape_spin(3),shape_spin(3))+1
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
        map_int(:,i_x,i_y)=map_int(:,i_x,i_y)+(sorien(spins(X:Z,i_x,i_y,i_m),spins(X:Z,ipu,i_y,i_m),spins(X:Z,ipu,ipv,i_m))+ &
         sorien(spins(X:Z,i_x,i_y,i_m),spins(X:Z,ipu,ipv,i_m),spins(X:Z,i_x,ipv,i_m)))/sqrt(2.0d0)

!!!!!!!!!!!!!!!!!!
! bottom surface
!!!!!!!!!!!!!!!!!!
        i_m=1

        map_int(:,i_x,i_y)=map_int(:,i_x,i_y)+(sorien(spins(X:Z,i_x,i_y,i_m),spins(X:Z,ipu,i_y,i_m),spins(X:Z,ipu,ipv,i_m))+ &
         sorien(spins(X:Z,i_x,i_y,i_m),spins(X:Z,ipu,ipv,i_m),spins(X:Z,i_x,ipv,i_m)))/sqrt(2.0d0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! unitl the end, i_m=1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        i_m=1
!!!!!!!!!!!!!!!!!!
! (0-10) surface (make a drawing...) the rotation is from u to w
!!!!!!!!!!!!!!!!!!
        ipm=2

        map_int(:,i_x,i_y)=map_int(:,i_x,i_y)+(sorien(spins(X:Z,i_x,i_y,i_m),spins(X:Z,ipu,i_y,i_m),spins(X:Z,ipu,i_y,ipm))+ &
         sorien(spins(X:Z,i_x,i_y,i_m),spins(X:Z,ipu,i_y,ipm),spins(X:Z,i_x,i_y,ipm)))/sqrt(2.0d0)
!!!!!!!!!!!!!!!!!!
! (100) surface (make a drawing...) the rotation is from v to w
!!!!!!!!!!!!!!!!!!
        map_int(:,i_x,i_y)=map_int(:,i_x,i_y)+(sorien(spins(X:Z,ipu,i_y,i_m),spins(X:Z,ipu,ipv,i_m),spins(X:Z,ipu,ipv,ipm))+ &
         sorien(spins(X:Z,ipu,i_y,i_m),spins(X:Z,ipu,ipv,ipm),spins(X:Z,ipu,i_y,ipm)))/sqrt(2.0d0)
!!!!!!!!!!!!!!!!!!
! (010) surface (make a drawing...) the rotation is from -u to w
!!!!!!!!!!!!!!!!!!
        map_int(:,i_x,i_y)=map_int(:,i_x,i_y)+(sorien(spins(X:Z,ipu,ipv,i_m),spins(X:Z,i_x,ipv,i_m),spins(X:Z,i_x,ipv,ipm))+ &
         sorien(spins(X:Z,ipu,ipv,i_m),spins(X:Z,i_x,ipv,ipm),spins(X:Z,ipu,ipv,ipm)))/sqrt(2.0d0)
!!!!!!!!!!!!!!!!!!
! (-100) surface (make a drawing...) the rotation is from w to v
!!!!!!!!!!!!!!!!!!
        map_int(:,i_x,i_y)=map_int(:,i_x,i_y)+(sorien(spins(X:Z,i_x,i_y,i_m),spins(X:Z,i_x,i_y,ipm),spins(X:Z,i_x,ipv,ipm))+ &
         sorien(spins(X:Z,i_x,i_y,i_m),spins(X:Z,i_x,ipv,ipm),spins(X:Z,i_x,ipv,i_m)))/sqrt(2.0d0)

       enddo
      enddo


#ifdef CPP_OPENMP
!$OMP end parallel do
#endif

#ifdef CPP_MPI
      if (i_ghost) then
       map_qorien=allreduce(map_int,shape_spin(2:3),3,MPI_COMM_BOX)
       else
       map_qorien=map_int
      endif
#endif

      map_qorien=map_int

      end subroutine qorien_2D_SL

      end module m_qorien
