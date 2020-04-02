module m_fftw

!
! kmesh
!
integer :: N_kpoint(3)
real(kind=8), allocatable :: kmesh(:,:)

!
! FFT of the dipolar martix D(r)
!
complex(kind=16), allocatable :: FFT_pos_D(:,:,:)

interface calculate_fft
  module procedure calculate_fft_vec_point,calculate_FFT_matrix
end interface

interface get_FFT
  module procedure get_FFT_vec_point,get_FFT_matrix
end interface

interface plot
  module procedure plot_2D,plot_3D
end interface

private
public :: get_k_mesh,calculate_fft
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

subroutine calculate_fft_vec_point(field,sense,dim_mode,astar,N)
use m_derived_types, only : vec_point
use m_get_position
implicit none
type(vec_point), intent(in) :: field(:)
real(kind=8), intent(in) :: sense       ! should be + or - one
real(kind=8), intent(in) :: astar(3,3)
integer, intent(in) :: dim_mode
! internal
integer :: i,Nsize,j,k,io_out,N_k,N(3)
complex(kind=16), allocatable :: FFT(:,:)
real(kind=8), allocatable :: pos(:,:)
real(kind=8) :: kvec(3)

kvec=0.0d0

Nsize=size(field)
N_k=size(kmesh,2)
allocate(FFT(dim_mode,N_k),pos(3,Nsize))
FFT=0.0d0
pos=0.0d0
call get_position(pos,'positions.dat')

do i=1,N_k
  kvec=matmul(kmesh(:,i),astar)

  FFT(:,i)=get_FFT(kvec,sense,pos,field,Nsize,dim_mode,N)/real(N_k)

enddo

call plot('FFT.dat',FFT)

end subroutine

!
! calculate the Fourrier transform of a matrix of reals
!

subroutine calculate_fft_matrix(field,sense,dim_lat,astar)
use m_vector, only : norm
implicit none
real(kind=8), intent(inout) :: field(:,:)
real(kind=8), intent(in) :: sense       ! should be + or - one
real(kind=8), intent(in) :: astar(:,:)
integer, intent(in) :: dim_lat(:)
! internal
integer :: Nx,Ny,Nz,Nksize,i,j,k,l,Nkx,Nky,Nkz
complex(kind=16), allocatable :: FFT(:,:,:,:,:)
real(kind=8), allocatable :: pos(:,:,:,:,:),real_dist(:,:,:,:)
real(kind=8) :: r(3),kvec(3)

Nx=dim_lat(1)
Ny=dim_lat(2)
Nz=dim_lat(3)

Nksize=product(N_kpoint)
Nkx=N_kpoint(1)
Nky=N_kpoint(2)
Nkz=N_kpoint(3)
allocate(FFT(3,3,Nkx,Nky,Nkz),pos(3,3,Nx,Ny,Nz),FFT_pos_D(3,3,Nksize),real_dist(3,Nx,Ny,Nz))
FFT=0.0d0
pos=0.0d0
FFT_pos_D=0.0d0
real_dist=0.0d0

real_dist=reshape(field,(/3,Nx,Ny,Nz/))

do i=1,Nz
  do j=1,Ny
    do k=1,Nx
      r=real_dist(:,k,j,i)
      pos(:,1,k,j,i)=(/ r(1)**2 - sum(r**2)/3 , r(1)*r(2) , r(1)*r(3) /)
      pos(:,2,k,j,i)=(/ r(1)*r(2) , r(2)**2 - sum(r**2)/3 , r(2)*r(3) /)
      pos(:,3,k,j,i)=(/ r(1)*r(3) , r(2)*r(3) , r(3)**2 - sum(r**2)/3 /)
    enddo
  enddo
enddo

l=0
do i=1,N_kpoint(3)
  do j=1,N_kpoint(2)
    do k=1,N_kpoint(1)
      l=l+1
      kvec=matmul(kmesh(:,l),astar)
      FFT(:,:,k,j,i)=get_FFT(pos,sense,kvec,real_dist)
    enddo
  enddo
enddo

FFT_pos_D=reshape(FFT,(/3,3,Nksize/))

call plot('FFT.dat',FFT_pos_D)

stop

end subroutine








!
! function that calculate the FFT of a matrix of vec_point for one component
!

function get_FFT_vec_point(kvec,sense,pos,field,Nsize,dim_mode,N)
use m_derived_types, only : vec_point
implicit none
integer, intent(in) :: Nsize,dim_mode,N(:)
real(kind=8), intent(in) :: kvec(:),pos(:,:),sense
type(vec_point), intent(in) :: field(:)
complex(kind=16) :: get_FFT_vec_point(dim_mode)
! internal
real(kind=8) :: alpha,r(3),Nx,Ny,Nz
integer :: k,j

get_FFT_vec_point=0.0d0
Nx=real(N(1))
Ny=real(N(2))
Nz=real(N(3))

do j=1,Nsize
    r(1)=pos(1,j)/Nx
    r(2)=pos(2,j)/Ny
    r(3)=pos(3,j)/Nz
    alpha=dot_product(kvec,r)

  do k=1,dim_mode

    get_FFT_vec_point(k)=get_FFT_vec_point(k)+complex(field(j)%w(k)*cos(sense*alpha),field(j)%w(k)*sin(sense*alpha))

  enddo
enddo

end function get_FFT_vec_point

!
! function that calculate the FFT of a matrix of real for one component
!

function get_FFT_matrix(pos,sense,kvec,real_dist)
implicit none
real(kind=8), intent(in) :: pos(:,:,:,:,:),sense,kvec(:),real_dist(:,:,:,:)
complex(kind=16) :: get_FFT_matrix(3,3)
! internal
integer :: shape_pos(5),i,j,k,Nx,Ny,Nz,l,m,Nsize
real(kind=8) :: alpha,r(3)

shape_pos=shape(pos)
Nx=shape_pos(3)
Ny=shape_pos(4)
Nz=shape_pos(5)
get_FFT_matrix=0.0d0
Nsize=product(shape_pos(3:5))

do i=1,Nz
  do j=1,Ny
    do k=1,Nx

      r=real_dist(:,k,j,i)
      alpha=dot_product(kvec,r)

        do l=1,3
          do m=1,3
            get_FFT_matrix(m,l)=get_FFT_matrix(m,l)+complex(pos(m,l,k,j,i)*cos(sense*alpha),pos(m,l,k,j,i)*sin(sense*alpha))
          enddo
        enddo

    enddo
  enddo
enddo

get_FFT_matrix=get_FFT_matrix/real(Nsize)

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
integer :: io_input
integer :: Nkpoint,i,iomp
logical :: i_plot,i_file

N_kpoint=my_lattice%dim_lat
i_plot=.false.
io_input=open_file_read(fname)
call get_parameter(io_input,fname,'kmesh',3,N_kpoint)
call get_parameter(io_input,fname,'write_kmesh',i_plot)
call close_file(fname,io_input)

Nkpoint=product(N_kpoint)
allocate(kmesh(3,Nkpoint))
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

end subroutine get_k_mesh




!
! plot the FFT
!

subroutine plot_2D(fname,FFT)
use m_io_utils
use m_io_files_utils
use m_convert
implicit none
character(len=*), intent(in) :: fname
complex(kind=16), intent(in) :: FFT(:,:)
! internal
integer :: io_out,shape_FFT(2),i,k
character(len=50) :: form

shape_FFT=shape(FFT)
form=convert('(',2*shape_FFT(1),'(f16.8,2x))')
io_out=open_file_write(fname)

do i=1,shape_FFT(2)
  write(io_out,form) (FFT(k,i),k=1,shape_FFT(1))
enddo

call close_file(fname,io_out)
end subroutine plot_2D

subroutine plot_3D(fname,FFT)
use m_io_utils
use m_io_files_utils
use m_convert
implicit none
character(len=*), intent(in) :: fname
complex(kind=16), intent(in) :: FFT(:,:,:)
! internal
integer :: io_out,shape_FFT(3),i,k,l
character(len=50) :: form

shape_FFT=shape(FFT)
form=convert('(',2*shape_FFT(1)*shape_FFT(2),'(f16.8,2x))')
io_out=open_file_write(fname)

do i=1,shape_FFT(3)
  write(io_out,form) ((FFT(k,l,i),k=1,shape_FFT(1)),l=1,shape_FFT(2))
enddo

call close_file(fname,io_out)
end subroutine plot_3D

end module m_fftw
