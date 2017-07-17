       module m_SkX_utils

       contains

!!! find the positions of all the skyrmions in the lattice

       subroutine find_XYsky(XSky,YSky,Nadd,qskx,dim_lat,net)
       use m_vector, only : norm
       implicit none
       integer, intent(in) :: Nadd,dim_lat(:)
       real(kind=8), intent(in) :: qskx,net(:,:)
       real(kind=8), intent(inout) :: XSky(Nadd),YSky(Nadd)
! internal
       integer :: i_x,i_y,N_x,N_y
       integer :: k

       N_x=nint(qskx*dim_lat(1))
       N_y=nint(qskx*dim_lat(2))

       k=1
       do i_x=1,N_x
        do i_y=1,N_y
        XSky(k)=dble(dim_lat(1)*(i_x-1))/dble(N_x)*net(1,1)+dble(dim_lat(2)*(i_y-1))/dble(N_y)*net(2,1)
        YSky(k)=dble(dim_lat(1)*(i_x-1))/dble(N_x)*net(1,2)+dble(dim_lat(2)*(i_y-1))/dble(N_y)*net(2,2)
         k=k+1
        enddo
       enddo

       if (k.ne.Nadd+1) then
        write(6,'(a)') 'problem in the lattice of skyrmion'
        write(6,'(a)') 'check SkX_utils routine'
        stop
       endif

       end subroutine find_XYsky

       end module m_SkX_utils
