      module m_fft
       interface fft
        module procedure fft_standard
#ifndef CPP_BRUTDIP
        module procedure fft_depol_r2c
        module procedure fft_depol_c2r
#endif
       end interface fft
      contains

!!!!!!!!!!!!!!!!!!!!!!
! easy FFT
!!!!!!!!!!!!!!!!!!!!!!
      subroutine fft_standard(dim_lat,kt)
      use m_rw_lattice, only : net
      use m_lattice, only : spin
      use m_vector, only : cross,norm
      use m_constants
!#ifdef CPP_OPENMP
!      use m_omp
!#endif
      implicit none
      integer, intent(in) :: dim_lat(3)
      real(kind=8), intent(in) :: kt
!internal
      double precision :: rx,ry
      complex*16 :: cdum1
      complex*16, allocatable :: fftcoef(:,:,:,:)
      integer :: i,j,i_lat,j_lat,k,i1,i2,j1,j2,qnx,qny,fin,i3,N_site
      real(kind=8) :: dum, kv(3,3),fft_norm(3,dim_lat(3)+1),kv0(3,3)
      character(len=30) :: fname,toto
      logical :: exists

      N_site=dim_lat(1)*dim_lat(2)*dim_lat(3)

      inquire (file='inp',exist=exists)
      if (.not. exists) then
      write(6,*) 'qnx and qny should be in inp'
      STOP
      endif

      open (121,file='inp',form='formatted',status='old',action='read')

      rewind(121)
      do
      read (121,'(a)',iostat=fin) toto
        if (fin /= 0) exit
        toto= trim(adjustl(toto))
        if (len_trim(toto)==0) cycle
        if ( toto(1:7) == 'gra_fft') then
           backspace(121)
           read(121,*) fname, exists, qnx, qny
          endif
      enddo
      close(121)

      allocate(fftcoef(qnx,qny,3,dim_lat(3)+1))
      fftcoef=dcmplx(0.d0,0.d0)
      fft_norm=0.0d0

      kv0(1,:) = pi(2.0d0)*cross(net(2,:),net(3,:))/dot_product(net(1,:),cross(net(2,:),net(3,:)))
      kv0(2,:) = pi(2.0d0)*cross(net(3,:),net(1,:))/dot_product(net(1,:),cross(net(2,:),net(3,:)))
      kv0(3,:) = pi(2.0d0)*cross(net(1,:),net(2,:))/dot_product(net(1,:),cross(net(2,:),net(3,:)))

      do i3=1,dim_lat(3)
       do i=1,dim_lat(1)
        do j=1,dim_lat(2)
         fft_norm(:,i3)=fft_norm(:,i3)+sum(Spin(4:6,i,j,i3,:))**2
        enddo
       enddo
       fft_norm(:,i3)=dsqrt(fft_norm(:,i3)/dble(dim_lat(1)*dim_lat(2)))
      enddo

#ifdef CPP_OPENMP
!$OMP parallel do default (shared) schedule(static) PRIVATE(i1,i2,j1,j2)
#endif

      do i3= 1,dim_lat(3)
      do i1= 1, qnx
       do i2= 1, qny

        kv(1,:) = kv0(1,:)*dble(i1-1)/dble(qnx-1)/2.0d0
        kv(2,:) = kv0(2,:)*dble(i2-1)/dble(qny-1)/2.0d0
        
        do j1 = 1, dim_lat(1)
         do j2 = 1, dim_lat(2)

         dum = dot_product(Spin(1:3,j1,j2,i3,1),kv(1,:)+kv(2,:))
         cdum1 = dcmplx(dcos(dum),dsin(dum))

         fftcoef(i1,i2,:,i3)=fftcoef(i1,i2,:,i3)+sum(Spin(4:6,j1,j2,i3,:))*cdum1
         enddo
        enddo

        fftcoef(i1,i2,:,i3) = fftcoef(i1,i2,:,i3)/dble(N_site)/fft_norm(:,i3)
       enddo
      enddo
      enddo

#ifdef CPP_OPENMP
!$OMP end parallel do
#endif

      do i3= 1,dim_lat(3)
      write(fname,'(f8.4,a,i1)') kT/k_B,'_',i3
      toto=trim(adjustl(fname))
      write(fname,'(a,18a,a)')'fft_',(toto(i:i),i=1,len_trim(toto)),'.dat'

      open(10+i3,file=fname,form='formatted')
      do j_lat=1,qnx
        do i_lat=1,qny
        write(10+i3,'(5(2x,f20.15))') norm(kv0(1,:)*dble(j_lat-1)/dble(qnx-1)/2.0d0), &
        norm(kv0(2,:)*dble(i_lat-1)/dble(qny-1)/2.0d0),&
        (dble(fftcoef(j_lat,i_lat,k,i3)**2+aimag(fftcoef(j_lat,i_lat,k,i3))**2),k=1,3)
        enddo
      if (mod(j_lat,qnx).eq.0) then
        write(10+i3,*) ''
      endif
      enddo
      close(10+i3)
      enddo

      deallocate(fftcoef)

      end subroutine fft_standard

#ifndef CPP_BRUTDIP
!!!!!!!!!!!!!!!!!!!!!!
! FFT dipole
!!!!!!!!!!!!!!!!!!!!!!
      subroutine fft_depol_r2c(Nx,Ny,Nz,rmat,cmat,rtrans,ctrans,plan)
      use m_rw_lattice, only : dim_lat
      use m_vector, only : cross,norm
      use, intrinsic :: iso_c_binding
      implicit none
      integer, intent(in) :: Nx,Ny,Nz
      real(kind=8), intent(in) :: rmat(:,:,:,:)
      complex, intent(out) :: cmat(:,:,:,:)
      type(C_PTR), intent(in) :: plan
      real(c_double), intent(inout) :: rtrans(:,:,:)
      complex(c_double_complex), intent(inout) :: ctrans(:,:,:)
!internal
      integer :: i,j,k,l
      integer :: N

      include 'fftw3.f03'

      N=size(rmat,1)

      ctrans=dcmplx(0.d0,0.d0)
      rtrans=0.0d0

      do i=1,N

       rtrans(1:Nx,1:Ny,1:Nz)=rmat(i,1:Nx,1:Ny,1:Nz)

       call fftw_execute_dft_r2c(plan,rtrans,ctrans)

       cmat(i,1:Nx,1:Ny,1:Nz)=ctrans(1:Nx,1:Ny,1:Nz)

      enddo

      end subroutine fft_depol_r2c

!!!!!!!!!!!!!!!!!!!!!!
! FFT dipole
!!!!!!!!!!!!!!!!!!!!!!
      subroutine fft_depol_c2r(Nx,Ny,Nz,cmat,rmat,alpha,rtrans,ctrans,plan)
      use m_rw_lattice, only : net,dim_lat
      use m_vector, only : cross,norm
      use, intrinsic :: iso_c_binding
      implicit none
      integer, intent(in) :: Nx,Ny,Nz
      integer,intent(in) :: alpha
      real(kind=8), intent(out) :: rmat(:,:,:,:)
      complex, intent(in) :: cmat(:,:,:,:)
      type(C_PTR), intent(in) :: plan
      real(c_double), intent(inout) :: rtrans(:,:,:)
      complex(c_double_complex), intent(inout) :: ctrans(:,:,:)
!internal
      integer :: i,j,k,l
      integer :: N

      include 'fftw3.f03'

      ctrans=dcmplx(0.d0,0.d0)
      rtrans=0.0d0

      N=size(cmat,1)

      do i=1,N

       ctrans(1:Nx,1:Ny,1:Nz)=cmat(i,1:Nx,1:Ny,1:Nz)

       call fftw_execute_dft_c2r(plan,ctrans,rtrans)

       rmat(i,1:Nx,1:Ny,1:Nz)=rtrans(1:Nx,1:Ny,1:Nz)/dble(Nx*Ny*Nz)
      enddo

!      call fftw_destroy_plan(plan)

      end subroutine fft_depol_c2r
#endif
      end module m_fft
