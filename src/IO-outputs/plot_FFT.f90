module m_plot_FFT
!HAS TO BE UPDATED TO NEW ORDER PARAMETER
#if 0

interface plot
  module procedure plot_2D,plot_3D
end interface

private
public :: plot_fft

contains

!
! calculate the Fourrier transform of a field of pointers for a particular iteration
!

!subroutine calculate_fft_vec_point_tag(field,sense,dim_mode,tag)
!use m_derived_types, only : vec_point
!use m_get_position
!use m_convert
!implicit none
!type(vec_point), intent(in) :: field(:)
!real(kind=8), intent(in) :: sense       ! should be + or - one
!integer, intent(in) :: dim_mode,tag
!! internal
!integer :: i,Nsize,N_k
!complex(kind=8), allocatable :: FFT(:,:)
!character(len=50) :: fname
!
!fname=convert('FFT_',tag,'.dat')
!Nsize=size(field)
!N_k=size(kmesh,2)
!allocate(FFT(dim_mode,N_k))
!FFT=0.0d0
!
!do i=1,N_k
!
!  FFT(:,i)=get_FFT(kmesh(:,i),sense,pos,field,Nsize,dim_mode)
!
!enddo
!
!call plot(fname,FFT)
!
!end subroutine

subroutine plot_fft(all_mode,sense,r,dim_lat,boundary,dim_mode,tag)
use m_derived_types, only : vec_point
use m_get_position
use m_fftw
use m_convert
implicit none
type(vec_point), intent(in) :: all_mode(:)
integer, intent(in) :: tag,dim_lat(:),dim_mode
real(kind=8), intent(in) :: sense,r(:,:)
logical, intent(in) :: boundary(:)
! internal
integer :: N
real(kind=8), allocatable :: positions(:,:),distances(:,:)
complex(kind=8), allocatable :: FFT(:,:)
character(len=50) :: fname

N=size(all_mode)
allocate(distances(3,N),positions(3,N),FFT(dim_mode,N))
distances=0.0d0
positions=0.0d0
FFT=0.0d0

call get_position(positions,'positions.dat')

!
! prepare the dipolar matrix (the matrix of the r to calculate the 1/r)
!

call calculate_distances(distances,positions,r,dim_lat,boundary)

call calculate_fft(all_mode,distances,sense,dim_mode,FFT)

fname=convert('FFT_',tag,'.dat')
call plot(fname,FFT)

end subroutine






!
! plot the FFT
!

subroutine plot_2D(fname,FFT)
use m_io_utils
use m_io_files_utils
use m_convert
implicit none
character(len=*), intent(in) :: fname
complex(kind=8), intent(in) :: FFT(:,:)
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
complex(kind=8), intent(in) :: FFT(:,:,:)
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
#endif
end module m_plot_FFT
