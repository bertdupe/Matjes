module m_kmesh
use m_io_files_utils

private
public :: get_kmesh

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Define a Monkhorst,Pack kmesh for Fourier transform and transport
! from PRB, 13, 5188 (1976)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_kmesh(dim_q,kv0,net,kmesh)
implicit none
integer,intent(in) :: dim_q(:)
real(kind=8), intent(in) :: kv0(:,:),net(:,:)
real(kind=8), intent(out) :: kmesh(:,:,:)
! internal
integer :: i,j,qnx,qny,io
real(kind=8) :: alpha(3),k(3)

qnx=dim_q(1)
qny=dim_q(2)
k=0.0d0

do i=1,3
    alpha(i)=dot_product(kv0(i,:),net(i,:))
enddo

io=open_file_write('kpoints')
do j=1,qny
   do i=1,qnx

   kmesh(1,i,j)=dble(2*i-qnx-1)/dble(2*qnx)
   kmesh(2,i,j)=dble(2*j-qny-1)/dble(2*qny)

   write(io,'(3f15.8,2x)') kmesh(:,i,j),0.0

   enddo
enddo

call close_file('kpoints',io)

end subroutine get_kmesh

end module m_kmesh
