module m_fftw

!
! kmesh in internal units
!
integer, protected, public :: N_kpoint(3)
real(kind=8), protected, public, allocatable :: kmesh(:,:)

!
! FFT of the dipolar martix D(r)
!
complex(kind=16), allocatable, protected, public :: FFT_pos_D(:,:,:)

interface calculate_fft
  module procedure calculate_fft_vec_point,calculate_FFT_matrix
  ! module procedure calculate_FFT_Hamiltonian
end interface

interface get_FFT
  module procedure get_FFT_vec_point,get_FFT_matrix,get_FFT_Dr_matrix
end interface

private
public :: get_k_mesh,calculate_fft,calculate_FFT_Dr
contains

#ifdef CPP_FFTW
!!!!!!!!!!!!!!!!!!!!!!
! FFT dipole
!!!!!!!!!!!!!!!!!!!!!!
subroutine fft_depol_r2c(Nx,Ny,Nz,rmat,cmat,rtrans,ctrans,plan)

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

!
! calculate the Fourrier transform of a field of pointers
!

subroutine calculate_fft_vec_point(field,pos,sense,dim_mode,FFT)
use m_derived_types, only : vec_point
implicit none
complex(kind=16), intent(inout) :: FFT(:,:)
type(vec_point), intent(in) :: field(:)
real(kind=8), intent(in) :: sense       ! should be + or - one
real(kind=8), intent(in) :: pos(:,:)
integer, intent(in) :: dim_mode
!! internal
integer :: i,Nsize,N_k

Nsize=size(field)
N_k=size(kmesh,2)
FFT=0.0d0

do i=1,N_k

  FFT(:,i)=get_FFT(kmesh(:,i),sense,pos,field,Nsize,dim_mode)

enddo

end subroutine

!
! calculate the Fourrier transform of a field of reals
!

subroutine calculate_fft_matrix(field,pos,sense,FFT)
use m_get_position
implicit none
real(kind=8), intent(in) :: field(:,:),pos(:,:)
real(kind=8), intent(in) :: sense       ! should be + or - one
complex(kind=16), intent(inout) :: FFT(:,:)
! internal
integer :: i,N_k,N(2)

N=shape(field)
N_k=size(kmesh,2)

do i=1,N_k

  FFT(:,i)=get_FFT(field,sense,kmesh(:,i),pos,N(1))

enddo

end subroutine









!
! calculate the Fourrier transform of the D(r) matrix
!

subroutine calculate_FFT_Dr(pos,sense,dim_lat)
use m_vector, only : norm
implicit none
real(kind=8), intent(in) :: sense       ! should be + or - one
real(kind=8), intent(in) :: pos(:,:)
integer, intent(in) :: dim_lat(:)
! internal
integer :: Nksize,Nsize,i
real(kind=8) :: r(3),norm_int
real(kind=8), allocatable :: pos_D(:,:,:)

Nksize=product(N_kpoint)
Nsize=product(dim_lat)
allocate(FFT_pos_D(3,3,Nksize),pos_D(3,3,Nsize))
FFT_pos_D=0.0d0

do i=1,Nsize

   r=pos(:,i)
   norm_int=norm(r)

   if (norm_int.gt.1.0d-8) then
      pos_D(:,1,i)=(/ r(1)**2 - sum(r**2)/3 , r(1)*r(2) , r(1)*r(3) /)/norm_int**5
      pos_D(:,2,i)=(/ r(1)*r(2) , r(2)**2 - sum(r**2)/3 , r(2)*r(3) /)/norm_int**5
      pos_D(:,3,i)=(/ r(1)*r(3) , r(2)*r(3) , r(3)**2 - sum(r**2)/3 /)/norm_int**5
   endif

enddo

do i=1,Nksize
      FFT_pos_D(:,:,i)=get_FFT(pos_D,sense,kmesh(:,i),pos,Nsize)
enddo

end subroutine








!
! function that calculate the FFT of a matrix of vec_point for one component
!

function get_FFT_vec_point(kvec,sense,pos,field,Nsize,dim_mode)
use m_derived_types, only : vec_point
implicit none
integer, intent(in) :: Nsize,dim_mode
real(kind=8), intent(in) :: kvec(:),pos(:,:),sense
type(vec_point), intent(in) :: field(:)
complex(kind=16) :: get_FFT_vec_point(dim_mode)
! internal
real(kind=8) :: phase,r(3)
integer :: k,j

get_FFT_vec_point=0.0d0

do j=1,Nsize
    r=pos(:,j)
    phase=dot_product(kvec,r)

  do k=1,dim_mode

    get_FFT_vec_point(k)=get_FFT_vec_point(k)+complex(field(j)%w(k)*cos(sense*phase),field(j)%w(k)*sin(sense*phase))

  enddo

enddo

get_FFT_vec_point=get_FFT_vec_point/real(Nsize)

end function get_FFT_vec_point

!
! function that calculate the FFT of the D(r) matrix of real for one component
!

function get_FFT_Dr_matrix(pos_D,sense,kvec,real_dist,Nsize)
implicit none
integer, intent(in) :: Nsize
real(kind=8), intent(in) :: pos_D(:,:,:),sense,kvec(:),real_dist(:,:)
complex(kind=16) :: get_FFT_Dr_matrix(3,3)
! internal
integer :: i,l,m
real(kind=8) :: phase,r(3)

get_FFT_Dr_matrix=0.0d0

do i=1,Nsize

   r=real_dist(:,i)
   phase=dot_product(kvec,r)

   do l=1,3
      do m=1,3
         get_FFT_Dr_matrix(m,l)=get_FFT_Dr_matrix(m,l)+complex(pos_D(m,l,i)*cos(sense*phase),pos_D(m,l,i)*sin(sense*phase))
      enddo
   enddo

enddo

get_FFT_Dr_matrix=get_FFT_Dr_matrix/real(Nsize)

end function get_FFT_Dr_matrix

!
! function that calculate the FFT of a matrix of real for one component
!

function get_FFT_matrix(field,sense,kvec,real_dist,dim_mode)
implicit none
real(kind=8), intent(in) :: field(:,:),sense,kvec(:),real_dist(:,:)
integer, intent(in) :: dim_mode
complex(kind=16) :: get_FFT_matrix(dim_mode)
! internal
integer :: shape_pos(2),i,l
real(kind=8) :: alpha,r(3)

shape_pos=shape(field)

get_FFT_matrix=0.0d0

do i=1,shape_pos(2)

   r=real_dist(:,i)
   alpha=dot_product(kvec,r)

   do l=1,shape_pos(1)
      get_FFT_matrix(l)=get_FFT_matrix(l)+complex(field(l,i)*cos(sense*alpha),field(l,i)*sin(sense*alpha))
   enddo

enddo

get_FFT_matrix=get_FFT_matrix/real(shape_pos(2))

end function get_FFT_matrix






!
! get the mesh for the Fourrier transform
!

subroutine get_k_mesh(fname,my_lattice)
use m_kmesh
use m_derived_types, only : lattice
use m_io_utils
use m_io_files_utils
implicit none
character(len=*), intent(in) :: fname
type(lattice), intent(in) :: my_lattice
! internal stuff
integer :: io_input,test
integer :: Nkpoint,i,iomp
logical :: i_plot,i_file

N_kpoint=my_lattice%dim_lat
i_plot=.false.
io_input=open_file_read(fname)
call get_parameter(io_input,fname,'kmesh',3,N_kpoint)
call get_parameter(io_input,fname,'write_kmesh',i_plot)
call close_file(fname,io_input)

Nkpoint=product(N_kpoint)

!
! check if the variable kmesh is allocated
!

allocate(kmesh(3,Nkpoint),stat=test)
if (test.ne.0) return

kmesh=0.0d0

i_file=.false.
inquire(file='kpoints',exist=i_file)
if (i_file) then
  io_input=open_file_read('kpoints')
  do iomp=1,Nkpoint
    read(io_input,*) (kmesh(i,iomp),i=1,3)
  enddo
  call close_file('kpoints',io_input)
else
  call get_kmesh(N_kpoint,kmesh,i_plot)
endif

!
! convert the kmesh in internal units
!

do iomp=1,Nkpoint
    kmesh(:,iomp)=matmul(kmesh(:,iomp),my_lattice%astar)
enddo

end subroutine get_k_mesh

end module m_fftw
