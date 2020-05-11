       module m_setup_dipole
       use, intrinsic :: iso_c_binding
       complex, dimension(:,:,:,:), allocatable :: ntensor,hcomplex,mcomplex
       real(kind=8), dimension(:,:,:,:), allocatable :: mmatrix,hreal
       integer :: Nfftx,Nffty,Nfftz
       type(C_PTR) :: planrtoc,planctor
       real(c_double),dimension(:,:,:), allocatable :: rtrans
       complex(c_double_complex),dimension(:,:,:), allocatable :: ctrans
       contains

! make the all the matrix of interest as power of 2 to prepare the FFT

      integer function findp2(N)
      implicit none
      integer :: N,p
      !internal
      integer :: i

      p=1
      do while (2**p.le.N)
       if (mod(N,2**p).eq.1) then
        exit
       endif
       findp2=2**p
       if (findp2.eq.N) return
       p=p+1
      enddo

      i=2
      do while (i.lt.N)
       i=i*2
      enddo

      findp2=i

      end function findp2

       end module m_setup_dipole

       subroutine setup_dipole(dim_lat,Periodic_log,net,n_mag_motif,dim_world)
       use m_setup_dipole
       use m_nmat
       implicit none
       integer, intent(in) :: dim_lat(3),n_mag_motif,dim_world
       real(kind=8), intent(in) :: net(3,3)
       logical, intent(in) :: Periodic_log(3)
!internal
! check the allocation of memory
      integer :: alloc_check

       include 'fftw3.f03'

       ! check that all dimensions are power of 2

       if ((Periodic_log(1).and.(.not.Periodic_log(2)).and.(.not.Periodic_log(3))).or.(dim_world.eq.1)) then
       Nfftx=findp2(dim_lat(1))
       Nffty=findp2(2*dim_lat(2))
       Nfftz=findp2(2*dim_lat(3)*n_mag_motif)
! case of a 2D system
       elseif (((.not.Periodic_log(3)).and.Periodic_log(2).and.Periodic_log(1)).or.(dim_world.eq.2)) then
       Nfftx=findp2(dim_lat(1))
       Nffty=findp2(dim_lat(2))
       Nfftz=findp2(2*dim_lat(3)*n_mag_motif)
       elseif ((.not.Periodic_log(1)).and.(.not.Periodic_log(2)).and.(.not.Periodic_log(3))) then
       Nfftx=findp2(2*dim_lat(1))
       Nffty=findp2(2*dim_lat(2))
       Nfftz=findp2(2*dim_lat(3)*n_mag_motif)
       else
       Nfftx=findp2(dim_lat(1))
       Nffty=findp2(dim_lat(2))
       Nfftz=findp2(dim_lat(3)*n_mag_motif)
       endif

!! allocation of the fourrier transform of the demag tensor
       allocate(ntensor(6,Nfftx,Nffty,Nfftz),stat=alloc_check)
       if (alloc_check.ne.0) write(6,'(a)') 'cannot allocate ntensor'
!! allocation of the real space magnetization matrix for the fft
       allocate(mmatrix(3,Nfftx,Nffty,Nfftz),stat=alloc_check)
       if (alloc_check.ne.0) write(6,'(a)') 'cannot allocate magnetization matrix'
!! allocation of the reciprocal space magnetization matrix for the fft
       allocate(mcomplex(3,Nfftx,Nffty,Nfftz),stat=alloc_check)
       if (alloc_check.ne.0) write(6,'(a)') 'cannot allocate reciprocal mag matrix'
!! allocation of the reciprocal space magnetization matrix for the fft
       allocate(hreal(3,Nfftx,Nffty,Nfftz),stat=alloc_check)
       if (alloc_check.ne.0) write(6,'(a)') 'cannot allocate real H matrix'
!! allocation of the reciprocal space magnetization matrix for the fft
       allocate(hcomplex(3,Nfftx,Nffty,Nfftz),stat=alloc_check)
       if (alloc_check.ne.0) write(6,'(a)') 'cannot allocate complex H matrix'
!! allocation of the real transfer matrix for the fft
       allocate(rtrans(Nfftx,Nffty,Nfftz),stat=alloc_check)
       if (alloc_check.ne.0) write(6,'(a)') 'cannot allocate real transfer matrix'
!! allocation of the complex transfer matrix for the fft
       allocate(ctrans(Nfftx,Nffty,Nfftz),stat=alloc_check)
       if (alloc_check.ne.0) write(6,'(a)') 'cannot allocate complex transfer matrix'

       mmatrix=0.0d0
       rtrans=0.0d0
       hreal=0.0d0
       ctrans=dcmplx(0.d0,0.d0)
       mcomplex=dcmplx(0.d0,0.d0)
       hcomplex=dcmplx(0.d0,0.d0)
       ntensor=dcmplx(0.d0,0.d0)

! prepare the plan for the FFT
       planrtoc = fftw_plan_dft_r2c_3d(Nfftz,Nffty,Nfftx,rtrans,ctrans,FFTW_ESTIMATE+FFTW_FORWARD)
       planctor = fftw_plan_dft_c2r_3d(Nfftz,Nffty,Nfftx,ctrans,rtrans,FFTW_ESTIMATE+FFTW_BACKWARD)

       if ((Periodic_log(1).and.(.not.Periodic_log(2)).and.(.not.Periodic_log(3))).or.(dim_world.eq.1)) then
        call Nmatrix(dim_lat,net(2,:),net(3,:),net,ntensor,Nfftx,Nffty,Nfftz,rtrans,ctrans,planrtoc)
       elseif (((.not.Periodic_log(3)).and.Periodic_log(2).and.Periodic_log(1)).or.(dim_world.eq.2)) then
        call Nmatrix(dim_lat,net(3,:),net,ntensor,Nfftx,Nffty,Nfftz,rtrans,ctrans,planrtoc)
       elseif ((.not.Periodic_log(1)).and.(.not.Periodic_log(2)).and.(.not.Periodic_log(3))) then
        call Nmatrix(dim_lat,net(1,:),net(2,:),net(3,:),net,ntensor,Nfftx,Nffty,Nfftz,rtrans,ctrans,planrtoc)
       else
        call Nmatrix(dim_lat,net,ntensor,Nfftx,Nffty,Nfftz,rtrans,ctrans,planrtoc)
       endif

       end subroutine setup_dipole
