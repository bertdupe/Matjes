! function that computes the demag matrix at each point of space for the convolution
! the arguments are
! u(3),v and w vectors along which you are not periodic. It must be a basis vector
! r(3,3) basis vector
! dim_lat
!
      module m_nmat
       interface Nmatrix
       module procedure und_nmat
       module procedure deuxd_nmat
       module procedure troid_nmat
       module procedure zerod_nmat
       end interface Nmatrix
      contains

! in case non periodicity in 1D only

      subroutine und_nmat(dim_lat,u,r,nmat,Nx,Ny,Nz,rtrans,ctrans,plan)
      use m_dten
      use m_fft
      use m_lattice, only : masque
#ifdef CPP_MPI
      use m_parameters, only : ierr,i_separate,i_average
      use m_mpi
      use m_mpi_prop, only : irank,isize,MPI_COMM
#endif
      use, intrinsic :: iso_c_binding
      implicit none
      integer, intent(in) :: dim_lat(3),Nx,Ny,Nz
      real(kind=8), intent(in) :: u(:),r(3,3)
      complex, intent(inout) :: nmat(:,:,:,:)
      type(C_PTR), intent(in) :: plan
      real(c_double), intent(inout) :: rtrans(:,:,:)
      complex(c_double_complex), intent(inout) :: ctrans(:,:,:)
      !dummy
      integer :: i,j,k,l,ii,jj,kk
      real(kind=8) :: vec(3),rij(3)
      real(kind=8) :: rtemp(6,Nx,Ny,Nz)
      integer :: N_start,N_stop

      rtemp=0.0d0

#ifdef CPP_MPI
      if ((.not.i_separate).and.(.not.i_average)) then
       N_start=dim_lat(2)/(isize)*irank+1
       N_stop=dim_lat(2)/(isize)*(irank+1)
       if (irank.eq.(isize-1)) N_stop=dim_lat(2)
       else
       N_start=1
       N_stop=dim_lat(2)
       endif
#else
       N_start=1
       N_stop=dim_lat(2)
#endif

      do k=1,dim_lat(3)
       do j=N_start,N_stop
        do i=1,dim_lat(1)

        if (masque(1,i,j,k).eq.0) cycle

         do kk=1-k,dim_lat(3)-k
          do jj=1,dim_lat(2)
           do ii=1,dim_lat(1)

          rij=dble(ii-1)*r(1,:)+dble(jj-1)*r(2,:)+ dble(kk)*r(3,:)

           rtemp(:,i,j,k)=rtemp(:,i,j,k)+demag_tensor(rij)

           enddo
          enddo
         enddo

        enddo
       enddo
      enddo

! non periodicity in 1D which is the z direction. Periodicity in x and y only
! the rtemp matrix has to be organised in bit reverse order. We start with z>dim_lat(3)
      do k=dim_lat(3)+2,2*dim_lat(3)
       do j=1,dim_lat(2)
        do i=1,dim_lat(1)
         rtemp(:,i,j,k)=rtemp(:,i,j,2*dim_lat(3)+2-k)
        enddo
       enddo
      enddo


#ifdef CPP_DEBUG
      do l=1,6
      do k=1,Nz
       write(*,*) 'Nz=',k
       do i=1,Nx
        write(*,*) (rtemp(l,i,j,k),j=1,Ny)
       enddo
      enddo
      enddo
#endif

#ifdef CPP_MPI
      rtemp=allreduce(rtemp,dim_lat,6,isize,MPI_COMM)
#endif

      call fft(Nx,Ny,Nz,rtemp,nmat,rtrans,ctrans,plan)

      end subroutine und_nmat

! in case non periodicity in 2D

      subroutine deuxd_nmat(dim_lat,u,v,r,nmat,Nx,Ny,Nz,rtrans,ctrans,plan)
      use m_dten
      use m_fft
      use m_lattice, only : masque
      use, intrinsic :: iso_c_binding
#ifdef CPP_MPI
      use m_parameters, only : ierr,i_separate,i_average
      use m_mpi
      use m_mpi_prop, only : irank,isize,MPI_COMM
#endif
      implicit none
      integer, intent(in) :: dim_lat(3),Nx,Ny,Nz
      real(kind=8), intent(in) :: u(:),v(:),r(3,3)
      complex, intent(inout) :: nmat(:,:,:,:)
      type(C_PTR), intent(in) :: plan
      real(c_double), intent(inout) :: rtrans(:,:,:)
      complex(c_double_complex), intent(inout) :: ctrans(:,:,:)
      !dummy
      real(kind=8) :: rtemp(6,Nx,Ny,Nz)
      integer :: i,j,k,ii,jj,kk,l
      real(kind=8) :: vec(3)
      integer :: N_start,N_stop

      rtemp=0.0d0
#ifdef CPP_MPI
      if ((.not.i_separate).and.(.not.i_average)) then
       N_start=dim_lat(2)/(isize)*irank+1
       N_stop=dim_lat(2)/(isize)*(irank+1)
       if (irank.eq.(isize-1)) N_stop=dim_lat(2)
       else
       N_start=1
       N_stop=dim_lat(2)
       endif
#else
       N_start=1
       N_stop=dim_lat(2)
#endif

! for each site
      do kk=1,dim_lat(3)
       do jj=N_start,N_stop
        do ii=1,dim_lat(1)

        if (masque(1,ii,jj,kk).eq.0) cycle

        do k=1-kk,dim_lat(3)-kk
         do j=1-jj,dim_lat(2)-jj
          do i=1,dim_lat(1)

         vec=dble(i-1)*r(1,:)+dble(j)*r(2,:)+dble(k)*r(3,:)

         rtemp(:,ii,jj,kk)=rtemp(:,ii,jj,kk)+demag_tensor(vec)

          enddo
         enddo
        enddo

        enddo
       enddo
      enddo

#ifdef CPP_MPI
      rtemp=allreduce(rtemp,dim_lat,6,isize,MPI_COMM)
#endif

! non periodicity in 2D which are the y and z direction. Periodicity in x only
! the rtemp matrix has to be organised in bit reverse order. We start with z<dim_lat(3)
      do k=1,dim_lat(3)

       do j=dim_lat(2)+2,2*dim_lat(2)
        do i=1,dim_lat(1)
         rtemp(:,i,j,k)=rtemp(:,i,2*dim_lat(2)+2-j,k)
        enddo
       enddo

      enddo

! the rtemp matrix has to be organised in bit reverse order. We start with z>dim_lat(3)
      do k=dim_lat(3)+2,2*dim_lat(3)

       do j=1,dim_lat(2)
        do i=1,dim_lat(1)
         rtemp(:,i,j,k)=rtemp(:,i,j,2*dim_lat(3)+2-k)
        enddo
       enddo

       do j=dim_lat(2)+2,2*dim_lat(2)
        do i=1,dim_lat(1)
         rtemp(:,i,j,k)=rtemp(:,i,2*dim_lat(2)+2-j,2*dim_lat(3)+2-k)
        enddo
       enddo

      enddo

#ifdef CPP_DEBUG
      do l=1,6
      do k=1,Nz
       do i=1,Nx
        write(6,*) (rtemp(l,i,j,k),j=1,Ny)
       enddo
       pause
      enddo
      enddo
#endif

      call fft(Nx,Ny,Nz,rtemp,nmat,rtrans,ctrans,plan)

      end subroutine deuxd_nmat

! in case non periodicity in 3D

      subroutine troid_nmat(dim_lat,u,v,w,r,nmat,Nx,Ny,Nz,rtrans,ctrans,plan)
      use m_dten
      use m_fft
      use m_lattice, only : masque
      use, intrinsic :: iso_c_binding
#ifdef CPP_MPI
      use m_parameters, only : ierr,i_separate,i_average
      use m_mpi
      use m_mpi_prop, only : irank,isize,MPI_COMM
#endif
      implicit none
      integer, intent(in) :: dim_lat(3),Nx,Ny,Nz
      real(kind=8), intent(in) :: u(:),v(:),w(:),r(3,3)
      complex, intent(inout) :: nmat(:,:,:,:)
      type(C_PTR), intent(in) :: plan
      real(c_double), intent(inout) :: rtrans(:,:,:)
      complex(c_double_complex), intent(inout) :: ctrans(:,:,:)
      !dummy
      real(kind=8) :: rtemp(6,Nx,Ny,Nz)
      integer :: i,j,k,ii,jj,kk,l
      real(kind=8) :: vec(3)
      integer :: N_start,N_stop

      rtemp=0.0d0

#ifdef CPP_MPI
      if ((.not.i_separate).and.(.not.i_average)) then
       N_start=dim_lat(2)/(isize)*irank+1
       N_stop=dim_lat(2)/(isize)*(irank+1)
       if (irank.eq.(isize-1)) N_stop=dim_lat(2)
       else
       N_start=1
       N_stop=dim_lat(2)
       endif
#else
       N_start=1
       N_stop=dim_lat(2)
#endif

      do k=1,dim_lat(3)
       do j=N_start,N_stop
        do i=1,dim_lat(1)

        if (masque(1,i,j,k).eq.0) cycle

        do kk=1-k,dim_lat(3)-k
         do jj=1-j,dim_lat(2)-j
          do ii=1-i,dim_lat(1)-i

         vec=dble(ii)*r(1,:)+dble(jj)*r(2,:)+dble(kk)*r(3,:)

         rtemp(:,i,j,k)=rtemp(:,i,j,k)+demag_tensor(vec)

          enddo
         enddo
        enddo

        enddo
       enddo
      enddo

#ifdef CPP_MPI
      rtemp=allreduce(rtemp,dim_lat,6,isize,MPI_COMM)
#endif
! the rtemp matrix has to be organised in bit reverse order. We start with z<dim_lat(3)
      do k=1,dim_lat(3)

       do j=1,dim_lat(2)
        do i=dim_lat(1)+2,2*dim_lat(1)
         rtemp(:,i,j,k)=rtemp(:,2*dim_lat(1)+2-i,j,k)
        enddo
       enddo

       do j=dim_lat(2)+2,2*dim_lat(2)
        do i=1,dim_lat(1)
         rtemp(:,i,j,k)=rtemp(:,i,2*dim_lat(2)+2-j,k)
        enddo
       enddo

       do j=dim_lat(2)+2,2*dim_lat(2)
        do i=dim_lat(2)+2,2*dim_lat(1)
         rtemp(:,i,j,k)=rtemp(:,2*dim_lat(1)+2-i,2*dim_lat(2)+2-j,k)
        enddo
       enddo

      enddo

! the rtemp matrix has to be organised in bit reverse order. We start with z>dim_lat(3)
      do k=dim_lat(3)+2,2*dim_lat(3)

       do j=1,dim_lat(2)
        do i=1,dim_lat(1)
         rtemp(:,i,j,k)=rtemp(:,i,j,2*dim_lat(3)+2-k)
        enddo
       enddo

       do j=1,dim_lat(2)
        do i=dim_lat(1)+2,2*dim_lat(1)
         rtemp(:,i,j,k)=rtemp(:,2*dim_lat(1)+2-i,j,2*dim_lat(3)+2-k)
        enddo
       enddo

       do j=dim_lat(2)+2,2*dim_lat(2)
        do i=1,dim_lat(1)
         rtemp(:,i,j,k)=rtemp(:,i,2*dim_lat(2)+2-j,2*dim_lat(3)+2-k)
        enddo
       enddo

       do j=dim_lat(2)+2,2*dim_lat(2)
        do i=dim_lat(2)+2,2*dim_lat(1)
         rtemp(:,i,j,k)=rtemp(:,2*dim_lat(1)+2-i,2*dim_lat(2)+2-j,2*dim_lat(3)+2-k)
        enddo
       enddo

      enddo

#ifdef CPP_DEBUG
      do l=1,6
      write(*,*) 'l=',l
      do k=1,Nz
       write(*,*) 'z=',k
       do i=1,Nx
        write(*,*) (rtemp(l,i,j,k),j=1,Ny)
       enddo
      enddo
      enddo
#endif

      call fft(Nx,Ny,Nz,rtemp,nmat,rtrans,ctrans,plan)

#ifdef CPP_DEBUG
      do l=1,6
      write(*,*) 'l=',l
      do k=1,Nz
       write(*,*) 'z=',k
       do i=1,Nx
        write(*,*) (abs(nmat(l,i,j,k)),j=1,Ny)
       enddo
      enddo
      enddo
#endif

      end subroutine troid_nmat

! in case of a periodic system

      subroutine zerod_nmat(dim_lat,r,nmat,Nx,Ny,Nz,rtrans,ctrans,plan)
      use m_dten
      use m_fft
      use m_lattice, only : masque
      use, intrinsic :: iso_c_binding
#ifdef CPP_MPI
      use m_parameters, only : ierr,i_separate,i_average
      use m_mpi
      use m_mpi_prop, only : irank,isize,MPI_COMM
#endif
      implicit none
      integer, intent(in) :: dim_lat(3),Nx,Ny,Nz
      real(kind=8), intent(in) :: r(3,3)
      complex, intent(inout) :: nmat(:,:,:,:)
      type(C_PTR), intent(in) :: plan
      real(c_double), intent(inout) :: rtrans(:,:,:)
      complex(c_double_complex), intent(inout) :: ctrans(:,:,:)
      !dummy
      real(kind=8) :: rtemp(6,Nx,Ny,Nz)
      integer :: i,j,k,ii,jj,kk,l
      real(kind=8) :: vec(3)
      integer :: N_start,N_stop

      rtemp=0.0d0

#ifdef CPP_MPI
      if ((.not.i_separate).and.(.not.i_average)) then
       N_start=dim_lat(2)/(isize)*irank+1
       N_stop=dim_lat(2)/(isize)*(irank+1)
       if (irank.eq.(isize-1)) N_stop=dim_lat(2)
       else
       N_start=1
       N_stop=dim_lat(2)
       endif
#else
       N_start=1
       N_stop=dim_lat(2)
#endif

      do k=1,dim_lat(3)
       do j=N_start,N_stop
        do i=1,dim_lat(1)

        if (masque(1,i,j,k).eq.0) cycle

        do kk=1,dim_lat(3)
         do jj=1,dim_lat(2)
          do ii=1,dim_lat(1)

           vec=dble(ii-1)*r(1,:)+dble(jj-1)*r(2,:)+dble(kk-1)*r(3,:)

           rtemp(:,i,j,k)=rtemp(:,i,j,k)+demag_tensor(vec)

          enddo
         enddo
        enddo

        enddo
       enddo
      enddo

#ifdef CPP_MPI
      rtemp=allreduce(rtemp,dim_lat,6,isize,MPI_COMM)
#endif

#ifdef CPP_DEBUG
      do l=1,6
      do k=1,Nz
       do i=1,Nx
        write(*,*) (rtemp(l,i,j,k),j=1,Ny)
       enddo
       pause
      enddo
      enddo
#endif

      call fft(Nx,Ny,Nz,rtemp,nmat,rtrans,ctrans,plan)

      end subroutine zerod_nmat

      end module m_nmat
