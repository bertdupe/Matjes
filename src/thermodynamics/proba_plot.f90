module m_proba_plot
use m_io_utils
use m_io_files_utils
use m_convert
use m_constants, only : k_B
use m_proba_base
implicit none

private
public :: histogram

contains

!!!!!!!!!!!!!!!!!
! print the probability distribution P(E)
!!!!!!!!!!!!!!!!!
subroutine histogram(this,tag)
class(proba_data), intent(in) :: this
integer, intent(in) :: tag

real(8) :: top,bottom,step
real(8), allocatable :: bins(:,:)
integer  :: io_out,i,N_bins
character(len=50) :: fname


N_bins=this%N_bins

allocate(bins(2,N_bins),source=0.0d0)
top = maxval(this%Pdistrib)
bottom = minval(this%Pdistrib)

step=(top-bottom)/real(N_bins)

do i=1,N_bins
   bins(1,i)=(bottom+step/2.0d0+step*real(i-1))/k_B
   bins(2,i)=count(((this%Pdistrib.gt.bottom+step*real(i-1)).and.(this%Pdistrib.le.bottom+step*real(i))))
enddo
bins(2,1)=bins(2,1)+1

fname=convert('histo_out_',tag)
io_out=open_file_write(fname)

do i=1,N_bins
   write(io_out,'(2(f12.6,2x))') bins(:,i)
enddo

call close_file(fname,io_out)

end subroutine

end module m_proba_plot
