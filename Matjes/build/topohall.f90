       module m_topohall
       interface topohall
        module procedure topohall_2D
       end interface topohall
       contains

! case of the calculation of the deviation angle depending of the position of the electron in the skyrmion
       subroutine topohall_2D(map_qorien,n,mag_motif)
       use m_vector, only : cross
       use m_constants, only : qel,pi
       use m_lattice, only : spin
       use m_derived_types, only : lattice
       implicit none
       real(kind=8), intent(in) :: map_qorien(:,:,:)
       type(lattice), intent(in) :: mag_motif
       integer, intent(in) :: n(:)

!internal
       real(kind=8) :: v(3),F(3),translation(3)
       real(kind=8) :: devia_F(2,n(1),n(1)),net(3,3)
       integer :: i_x,i_y,j

       v=0.0d0
       v(1)=1.0d0
       F=0.0d0
       translation=(net(1,:)+net(2,:)+net(3,:))/2.0d0
       net=mag_motif%areal

       devia_F=0.0d0

       do i_y=1,n(2)
        do i_x=1,n(1)
         F=cross(v,map_qorien(:,i_x,i_y))/pi(4.0d0)
        devia_F(1,i_x,i_y)=F(1)
        devia_F(2,i_x,i_y)=F(2)
        enddo
       enddo

! writting part
       write(6,'(a)') 'plot the topological charge'

       OPEN(70,FILE='force_Hall.dat',action='write',status='unknown',form='formatted')

       do i_y=1,n(2)
        do i_x=1,n(1)
         Write(70,'(5(E20.10E3,2x))') (Spin(j,i_x,i_y,1,1)+translation(j),j=1,3), &
     &    (devia_F(j,i_x,i_y),j=1,2)
        enddo
       enddo

       close(70)

       end subroutine topohall_2D

       end module m_topohall
