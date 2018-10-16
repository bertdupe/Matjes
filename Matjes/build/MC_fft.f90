       subroutine MC_fft(spins,pos,Im,Qim,dim_lat,N_cell,net)
       use m_lattice, only : masque
       use m_constants, only : pi
       use m_vector, only : cross,norm
#ifdef CPP_MPI
       use m_parameters, only : n_ghost
       use m_mpi_prop, only : irank_working
#endif
       implicit none
       integer, intent(in) :: dim_lat(3),N_cell
       real(kind=8), intent(in) :: spins(3,dim_lat(1),dim_lat(2),dim_lat(3)),pos(3,dim_lat(1),dim_lat(2),dim_lat(3))
       real(kind=8), intent(in) :: net(3,3)
       real(kind=8), intent(inout) :: Im(dim_lat(3)),Qim(dim_lat(3))
!!! fft part
      integer :: i1,i2,i3,j1,j2,qnx,qny,indmax(dim_lat(3),2)
      integer :: qnxmax(dim_lat(3)),qnymax(dim_lat(3)),qnzmax
      integer :: q1(dim_lat(3)),q2(dim_lat(3))
      real(kind=8) :: dum, kv(3,3), kv0(3,3),fftr_S,ffti_S,dum_norm
      real(kind=8), allocatable :: fftcoef(:,:,:)
      logical :: exists,exi_fft
!! FFT of the norm of the spins
      real(kind=8) :: norm_S(dim_lat(1),dim_lat(2),dim_lat(3))
! dummy
      integer :: ix,iy,iz,io
#ifdef CPP_MPI
      character(len=50) :: fname
      character(len=26) :: toto
      integer :: ierr(3)

      include 'mpif.h'
#endif

      q1=1
      q2=1
!! prepare the FFT
      kv0(1,:) = pi(2.0d0)*cross(net(2,:),net(3,:))/dot_product(net(1,:),cross(net(2,:),net(3,:)))
      kv0(2,:) = pi(2.0d0)*cross(net(3,:),net(1,:))/dot_product(net(1,:),cross(net(2,:),net(3,:)))
      kv0(3,:) = pi(2.0d0)*cross(net(1,:),net(2,:))/dot_product(net(1,:),cross(net(2,:),net(3,:)))
!!!!
#ifdef CPP_MPI
      write(fname,'(I6)') irank_working/n_ghost
      toto=trim(adjustl(fname))
      write(fname,'(a,14a,a)')'fft',(toto(i:i),i=1,len_trim(toto)),'.in'
      inquire (file=fname,exist=exists)
      call mpi_allreduce(exists,exi_fft,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,ierr)
#else
      inquire (file='fft.in',exist=exi_fft)
#endif
      exists=.False.
      if (exi_fft) then
#ifdef CPP_MPI
       io=120+irank_working
       open (io,file=fname,form='formatted',status='old', &
        action='read')
#else
       io=121
       open (io,file='fft.in',form='formatted',status='old', &
        action='read')
#endif
       rewind(io)
       do i3=1 ,dim_lat(3)
       read(io,*) qnxmax(i3), qnymax(i3)
       enddo
       close(io)
       q1=1
       q2=1
      else
       qnxmax=dim_lat(1)
       qnymax=dim_lat(2)
      endif

      inquire (file='Q_fft.in',exist=exists)
      if (exists) then
       open (121,file='Q_fft.in',form='formatted',status='old', &
        action='read')
       rewind(121)
       do i3=1 ,dim_lat(3)
       read(121,*) q1(i3), q2(i3)
       enddo
       close(121)
       qnxmax=q1+2
       qnymax=q2+2
      endif
      qnzmax=dim_lat(3)


      allocate(fftcoef(maxval(qnxmax,1),maxval(qnymax,1),qnzmax))


      fftcoef=0.d0
      dum_norm=0.0d0
      qnx=dim_lat(1)
      qny=dim_lat(2)

#ifdef CPP_OPENMP
!$OMP parallel DO REDUCTION(+:dum_norm) private(i) default(shared)
#endif
      do iz=1,dim_lat(3)
       do iy=1,dim_lat(2)
        do ix=1,dim_lat(1)
        if (masque(1,ix,iy,iz).eq.0) cycle
        norm_S(ix,iy,iz)=sqrt(spins(1,ix,iy,iz)**2+spins(2,ix,iy,iz)**2+spins(3,ix,iy,iz)**2)
        if (norm_S(ix,iy,iz).gt.1.0d-8) dum_norm=dum_norm+(spins(3,ix,iy,iz)/norm_S(ix,iy,iz))**2
        enddo
       enddo
      enddo
#ifdef CPP_OPENMP
!$OMP end parallel do
#endif

      do i3= 1, dim_lat(3)
      do i2= q2(i3), qnymax(i3)
       do i1= q1(i3), qnxmax(i3)

        kv(1,:) = kv0(1,:)*dble(i1-1)/dble(qnx-1)/2.0d0
        kv(2,:) = kv0(2,:)*dble(i2-1)/dble(qny-1)/2.0d0
        fftr_S=0.0d0
        ffti_S=0.0d0
#ifdef CPP_OPENMP
!$OMP parallel DO private(j1,j2,k) reduction(+:fftr_S,ffti_S) default(shared)
#endif
        do j2 = 1, dim_lat(2)
         do j1 = 1, dim_lat(1)
         if ((masque(1,j1,j2,i3).eq.0).or.(norm_S(j1,j2,i3).lt.1.0d-8)) cycle

         dum = dot_product(pos(:,j1,j2,i3),kv(1,:)+kv(2,:))

        fftr_S=fftr_S+spins(3,j1,j2,i3)/norm_S(j1,j2,i3)*dcos(dum)
        ffti_S=ffti_S+spins(3,j1,j2,i3)/norm_S(j1,j2,i3)*dsin(dum)

         enddo
        enddo
#ifdef CPP_OPENMP
!$OMP end parallel do
#endif
         fftcoef(i1,i2,i3) = dsqrt(fftr_S**2+ffti_S**2)/dum_norm/dble(N_cell)

       enddo
      enddo
      enddo

!! find the maximum in the coeff

      do i3= 1, dim_lat(3)
      Im(i3)=maxval(fftcoef(:,:,i3))
      indmax(i3,:)=maxloc(fftcoef(:,:,i3))
      Qim(i3)=norm(dble(indmax(i3,1)-1)*kv0(1,:)/dble(qnx-1)*norm(net(1,:))/pi(2.0d0)+ &
      dble(indmax(i3,2)-1)*kv0(2,:)/dble(qny-1))*norm(net(2,:))/pi(2.0d0)

      if (fftcoef(indmax(i3,1),indmax(i3,2),i3).lt.1.0d-10) indmax(i3,:)=(/1,1/)

      if ((indmax(i3,1).gt.dim_lat(1)/2).or.(indmax(i3,2).gt.dim_lat(2)/2)) &
       indmax(i3,:)=(/dim_lat(1),dim_lat(2)/)/2
      enddo

!! if fft.in is not set then set the window
      if (.not.(exi_fft).and.(.not.exists)) then
#ifdef CPP_MPI
       io=irank_working+105
       write(fname,'(I5)') irank_working
       toto=trim(adjustl(fname))
       write(fname,'(a,14a,a)')'fft',(toto(i:i),i=1,len_trim(toto)),'.in'
       open (io,file=fname,form='formatted',status='new',action='write')
       rewind(io)
       do i3= 1, qnzmax
       write(io,'(2(I8))') indmax(i3,1),indmax(i3,2)
       enddo
       close(io)
#else
       write(6,*) 'set fft.in'
       open (121,file='fft.in',form='formatted',status='unknown',action='write')
       rewind(121)
       do i3= 1, qnzmax
       write(121,'(2(I8))') indmax(i3,1),indmax(i3,2)
       enddo
       close(121)
#endif
      endif

      deallocate(fftcoef)
      end subroutine MC_fft
