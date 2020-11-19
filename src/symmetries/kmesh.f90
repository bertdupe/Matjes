module m_kmesh
use m_io_files_utils

private
public :: get_kmesh

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Define a Monkhorst,Pack kmesh for Fourier transform and transport
! from PRB, 13, 5188 (1976)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_kmesh(dim_q,kmesh,i_plot)
use m_get_position, only : get_position_ND_to_1D
implicit none
integer,intent(in) :: dim_q(:)
logical, intent(in) :: i_plot
real(kind=8), intent(inout) :: kmesh(:,:)
! internal
integer :: i,j,k,qnx,qny,qnz,io,iomp,Ilat(3)

qnx=dim_q(1)
qny=dim_q(2)
qnz=dim_q(3)

if (i_plot) io=open_file_write('kpoints')

do k=1,qnz
  do j=1,qny
    do i=1,qnx

      Ilat=(/i,j,k/)
      iomp=get_position_ND_to_1D(Ilat,dim_q)

      kmesh(1,iomp)=dble(2*i-qnx-1)/dble(2*qnx)
      kmesh(2,iomp)=dble(2*j-qny-1)/dble(2*qny)
      kmesh(3,iomp)=dble(2*k-qnz-1)/dble(2*qnz)

      if (i_plot) write(io,'(3f15.8,2x)') kmesh(:,iomp)

     enddo
  enddo
enddo

if (i_plot) call close_file('kpoints',io)

end subroutine get_kmesh

end module m_kmesh
