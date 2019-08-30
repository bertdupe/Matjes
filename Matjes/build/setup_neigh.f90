       module m_voisins
       interface periodic
        module procedure period
       end interface periodic
       contains

       subroutine period(ind,nvoi,d,NN,mask,tableNN,pos,my_lattice,my_motif)
       use m_vector , only : norm
       use m_derived_types, only : lattice,cell
#ifdef CPP_MPI
       use m_make_box, only : Xstart,Ystart,Zstart,Zstop,Ystop,Xstop
       use m_mpi_prop, only : MPI_COMM,irank,start
       use m_parameters, only : i_ghost
#endif
       implicit none
! variable that come in
       integer, intent(in) :: nvoi,NN,tableNN(:,:,:,:,:,:)
       integer, intent(in) :: ind(:)
       real(kind=8), intent(in) :: d(:),pos(:,:,:,:,:)
       type(lattice), intent(in) :: my_lattice
       type(cell), intent(in) :: my_motif
! value of the function
       integer, intent(inout) :: mask(:,:,:,:)
!dummy variable
       integer :: i_x,i_y,i_z,i_m,k,l,avant,i,dim_lat(3),nmag
       integer, allocatable :: masque(:,:,:,:)
       real (kind=8) :: dist,vec(3),net(3,3)
       logical :: exists,boundary(3)
! position of the neighbors
       integer :: v_x,v_y,v_z,v_m,ig_x,ig_y,ig_z
#ifndef CPP_MPI
       integer, parameter ::  Xstart=1
       integer, parameter ::  Ystart=1
       integer, parameter ::  Zstart=1
       integer ::  Xstop,Ystop,Zstop
       integer, dimension(3) :: start=0
#else
       integer :: n_comm,ierr(3)
#endif

#ifdef CPP_MPI
       include 'mpif.h'
#endif

       dim_lat=my_lattice%dim_lat
       net=my_lattice%areal
       nmag=size(my_motif%i_mom)
       boundary=my_lattice%boundary
       allocate(masque(nvoi,dim_lat(1),dim_lat(2),dim_lat(3)))

#ifdef CPP_MPI
       n_comm=product(dim_lat)*nvoi
#else
       Xstop=dim_lat(1)
       Ystop=dim_lat(2)
       Zstop=dim_lat(3)
#endif

       masque=0
       masque(1,Xstart:Xstop,Ystart:Ystop,Zstart:Zstop)=1
       avant=0
       exists=.False.

       inquire (file='wire.shape',exist=exists)
       if (exists) then

       do k=1,NN
#ifdef CPP_OPENMP
!$OMP parallel private(i_x,i_y,i_z,i_m,l,v_x,v_y,v_z,v_m,vec,dist) default(shared)
#endif
       do l=1,ind(k)
        do i_m=1,nmag
        if (.not.my_motif%i_mom(i_m)) cycle
         do i_z=Zstart,Zstop
          do i_y=Ystart,Ystop
           do i_x=Xstart,Xstop

           ig_x=i_x-start(1)
           ig_y=i_y-start(2)
           ig_z=i_z-start(3)

           v_x=tableNN(1,avant+l,ig_x,ig_y,ig_z,i_m)
           v_y=tableNN(2,avant+l,ig_x,ig_y,ig_z,i_m)
           v_z=tableNN(3,avant+l,ig_x,ig_y,ig_z,i_m)
           v_m=tableNN(4,avant+l,ig_x,ig_y,ig_z,i_m)

         vec=pos(1:3,i_x,i_y,i_z,i_m)-pos(1:3,v_x,v_y,v_z,v_m)

         dist=norm(vec)

          if ((abs(dist-d(k)).gt.1.0d-8).and.(dabs(dot_product(vec,net(2,:))).gt.d(k))) then
           masque(avant+l+1,i_x,i_y,i_z)=0
           else
           masque(avant+l+1,i_x,i_y,i_z)=1
          endif

           enddo
           enddo
          enddo
         enddo
        enddo
#ifdef CPP_OPENMP
!$OMP end parallel
#endif
        avant=avant+ind(k)
       enddo

! periodic boundary condition in a square cell. The motif is not cut
       else

       do k=1,NN
#ifdef CPP_OPENMP
!$OMP parallel private(i_x,i_y,i_z,i_m,l,v_x,v_y,v_z,v_m,vec,dist) default(shared)
#endif

       do l=1,ind(k)
        do i_m=1,nmag
        if (.not.my_motif%i_mom(i_m)) cycle
         do i_z=Zstart,Zstop
          do i_y=Ystart,Ystop
           do i_x=Xstart,Xstop

           ig_x=i_x-start(1)
           ig_y=i_y-start(2)
           ig_z=i_z-start(3)

           v_x=tableNN(1,avant+l,ig_x,ig_y,ig_z,i_m)
           v_y=tableNN(2,avant+l,ig_x,ig_y,ig_z,i_m)
           v_z=tableNN(3,avant+l,ig_x,ig_y,ig_z,i_m)
           v_m=tableNN(4,avant+l,ig_x,ig_y,ig_z,i_m)

       dist=norm(pos(1:3,i_x,i_y,i_z,i_m)-pos(1:3,v_x,v_y,v_z,v_m))
       if (abs(dist-d(k)).gt.1.0d-8) then
        do i=1,3
         if (boundary(i)) cycle
         dist=abs(pos(i,i_x,i_y,i_z,i_m)-pos(i,v_x,v_y,v_z,v_m))
         if (abs(dist-d(k)).gt.1.0d-8) masque(avant+l+1,i_x,i_y,i_z)=0
        enddo
       else
        masque(avant+l+1,i_x,i_y,i_z)=1
       endif

           enddo
           enddo
          enddo
         enddo
        enddo

#ifdef CPP_OPENMP
!$OMP end parallel
#endif
        avant=avant+ind(k)
       enddo

       endif

#ifdef CPP_MPI
      if (i_ghost)then
       call mpi_allreduce(masque,mask,N_comm,MPI_INTEGER,MPI_SUM,MPI_COMM,ierr)
       write(6,*) irank, " gather the masque"
      endif
#else
      mask=masque
#endif

#ifdef CPP_DEBUG
#ifdef CPP_MPI
      if (irank.eq.0) then
#endif
      do i_z=1,dim_lat(3)
       do i_y=1,dim_lat(2)
        do i_x=1,dim_lat(1)
        write(*,*) i_x,i_y,i_z,mask(:,i_x,i_y,i_z)
        enddo
       enddo
      enddo
      stop
#ifdef CPP_MPI
      endif
#endif
#endif

#ifdef CPP_MPI
       if (irank.eq.0) write(6,'(a)') "periodic boundary conditions setup is done"
#else
       write(6,'(a)') "periodic boundary conditions setup is done"
#endif
       end subroutine period

       end module
