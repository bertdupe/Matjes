module m_fft
use m_derived_types
use m_io_files_utils
use m_io_utils
use m_get_position
use m_kmesh

interface fft
    module procedure fft_modes
#ifndef CPP_BRUTDIP
    module procedure fft_depol_r2c
    module procedure fft_depol_c2r
#endif
end interface fft

private
public :: fft

contains

!!!!!!!!!!!!!!!!!!!!!!
! easy FFT
!!!!!!!!!!!!!!!!!!!!!!
subroutine fft_modes(my_lattice,my_motif,kt)
use m_vector, only : cross,norm
use m_constants

implicit none
type(lattice),intent(in) :: my_lattice
type(cell), intent(in) :: my_motif
real(kind=8),optional,intent(in) :: kt
!internal
complex(kind=16) :: cdum1
complex(kind=16), allocatable :: fftcoef(:,:,:,:)
real(kind=8), allocatable :: position(:,:,:,:,:)
integer :: i,j,i1,i2,j1,j2,qnx,qny,i3,N_site,dim_lat(3)
integer :: io
real(kind=8) :: dum, kv(3,3),kv0(3,3), net(3,3),kt_int
real(kind=8), allocatable :: fft_norm(:,:),kmesh(:,:,:)
character(len=30) :: fname,toto
logical :: exists

dim_lat=my_lattice%dim_lat
net=my_lattice%areal
kv0=my_lattice%astar
N_site=dim_lat(1)*dim_lat(2)*dim_lat(3)
qnx=dim_lat(1)
qny=dim_lat(2)

io=open_file_read('input')
call get_parameter(io,'input','gra_fft',exists)
! you do not want to do the FFT so leave the routine
if (.not.exists) then
  call close_file('input',io)
  return
endif

call get_parameter(io,'input','qnx',qnx)
call get_parameter(io,'input','qny',qny)

call close_file('input',io)

if (present(kt)) then
  kt_int=kt
else
  kt_int=0.0d0
endif

allocate(kmesh(2,qnx,qny))
call get_kmesh((/qnx,qny,1/),kv0,net,kmesh)

allocate(position(3,dim_lat(1),dim_lat(2),dim_lat(3),1))

call get_position(position,dim_lat,net,my_motif)

allocate(fft_norm(3,dim_lat(3)+1))

allocate(fftcoef(3,qnx,qny,dim_lat(3)+1))
fftcoef=dcmplx(0.d0,0.d0)
fft_norm=0.0d0

do i3=1,dim_lat(3)
   do j=1,dim_lat(2)
      do i=1,dim_lat(1)
        fft_norm(:,i3)=fft_norm(:,i3)+sum(my_lattice%l_modes(i,j,i3,1)%w(1:3)**2)
      enddo
   enddo
   fft_norm(:,i3)=sqrt(fft_norm(:,i3)/dble(dim_lat(1)*dim_lat(2)))
enddo

#ifdef CPP_OPENMP
!$OMP parallel do default (shared) schedule(static) PRIVATE(i1,i2,j1,j2)
#endif

do i3= 1,dim_lat(3)
   do i2= 1, qny
      do i1= 1, qnx

       kv(1,:) = kv0(1,:)*kmesh(1,i1,i2)
       kv(2,:) = kv0(2,:)*kmesh(2,i1,i2)
        
        do j2 = 1, dim_lat(2)
          do j1 = 1, dim_lat(1)

         dum = dot_product(position(:,j1,j2,i3,1),kv(1,:)+kv(2,:))
         cdum1 = dcmplx(dcos(dum),dsin(dum))

         fftcoef(:,i1,i2,i3)=fftcoef(:,i1,i2,i3)+my_lattice%l_modes(j1,j2,i3,1)%w(1:3)*cdum1
         enddo
        enddo

        fftcoef(:,i1,i2,i3) = fftcoef(:,i1,i2,i3)/fft_norm(:,i3)/dble(N_site)
      enddo
   enddo
enddo

#ifdef CPP_OPENMP
!$OMP end parallel do
#endif

do i3= 1,dim_lat(3)

  write(fname,'(f8.4,a,i1)') kt_int/k_B,'_',i3
  toto=trim(adjustl(fname))
  write(fname,'(a,18a,a)')'fft_',(toto(i:i),i=1,len_trim(toto)),'.dat'

  io=open_file_write(fname)
  call dump_config(io,kv0,kmesh,fftcoef(:,:,:,i3))

  call close_file(fname,io)
enddo

deallocate(fftcoef,position,kmesh)

end subroutine fft_modes

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
