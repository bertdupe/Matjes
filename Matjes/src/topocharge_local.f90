       module m_topocharge_local
       interface local_topo
        module procedure local_topo_2D
        module procedure local_topo_2D_SL
        module procedure local_topo_3D_SL
       end interface local_topo

       interface local_vortex
        module procedure local_vortex_2D
       end interface local_vortex

       contains

       subroutine local_topo_2D(i_x,i_y,qm,qp,spin,shape_spin,masque,shape_masque)
       use m_parameters, only : Periodic_log,i_print_W
       use m_vector, only : area
       implicit none
! input
      integer, intent(in) :: shape_spin(5),shape_masque(4)
      real(kind=8), intent(in) :: spin(shape_spin(1),shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5))
      integer, intent(in) :: masque(shape_masque(1),shape_masque(2),shape_masque(3),shape_masque(4))
       integer, intent(in) :: i_x,i_y
       real(kind=8), intent(out) :: qm,qp
! internal
       integer :: i,j,ipu,ipv,ii_x,ii_y
       real(kind=8) :: charge

       charge=0.0d0
       qm=0.0d0
       qp=0.0d0

       do j=-1,1,2
        do i=-1,1,2

        if (Periodic_log(1).and.Periodic_log(2)) then
         ipu=mod(i_x+i-1+shape_spin(2),shape_spin(2))+1
         ipv=mod(i_y+j-1+shape_spin(3),shape_spin(3))+1
        else
         if (i_x+i.eq.shape_spin(2)) then
          ipu=i_x
          else
          ipu=i_x+i
         endif
         if (i_y+j.eq.shape_spin(3)) then
          ipv=i_y
          else
          ipv=i_y+1
         endif
        endif

        charge=charge+(area(Spin(4:6,i_x,i_y,1,1),Spin(4:6,ipu,i_y,1,1),Spin(4:6,ipu,ipv,1,1))+ &
        area(Spin(4:6,i_x,i_y,1,1),Spin(4:6,ipu,ipv,1,1),Spin(4:6,i_x,ipv,1,1)))*dble(i*j)

        enddo
       enddo

       if (charge.gt.1.0d-10) qp=charge/4.0d0
       if (charge.lt.-1.0d-10) qm=charge/4.0d0

       end subroutine local_topo_2D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       subroutine local_topo_2D_SL(i_x,i_y,i_motif,qm,qp,spin,shape_spin,masque,shape_masque)
       use m_parameters, only : Periodic_log
       use m_vector, only : area
       implicit none
! input
       integer, intent(in) :: shape_spin(5),shape_masque(4)
       real(kind=8), intent(in) :: spin(shape_spin(1),shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5))
       integer, intent(in) :: masque(shape_masque(1),shape_masque(2),shape_masque(3),shape_masque(4))
       integer, intent(in) :: i_x,i_y,i_motif
       real(kind=8), intent(out) :: qm,qp
! internal
       integer :: i,j,ipu,ipv,i_m,ipm
       real(kind=8) :: charge,nmag_motif

       qm=0.0d0
       qp=0.0d0
       charge=0.0d0
       nmag_motif=real(shape_spin(5))

       do j=-1,1,2
        do i=-1,1,2

        if (Periodic_log(1).and.Periodic_log(2)) then
         ipu=mod(i_x+i-1+shape_spin(2),shape_spin(2))+1
         ipv=mod(i_y+j-1+shape_spin(3),shape_spin(3))+1
        else
         if (i_x+i.eq.shape_spin(2)) then
          ipu=i_x
          else
          ipu=i_x+i
         endif
         if (i_y+j.eq.shape_spin(3)) then
          ipv=i_y
          else
          ipv=i_y+1
         endif
        endif

         do i_m=1,shape_spin(5)

          charge=charge+(area(Spin(4:6,i_x,i_y,1,i_m),Spin(4:6,ipu,i_y,1,i_m),Spin(4:6,ipu,ipv,1,i_m))+ &
           area(Spin(4:6,i_x,i_y,1,i_m),Spin(4:6,ipu,ipv,1,i_m),Spin(4:6,i_x,ipv,1,i_m)))*dble(i*j)

         enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! unitl the end, i_m=1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         i_m=1
!!!!!!!!!!!!!!!!!!
! (0-10) surface (make a drawing...) the rotation is from u to w
!!!!!!!!!!!!!!!!!!
         ipm=2

         charge=charge+(area(Spin(4:6,i_x,i_y,1,i_m),Spin(4:6,ipu,i_y,1,i_m),Spin(4:6,ipu,i_y,1,ipm))+ &
          area(Spin(4:6,i_x,i_y,1,i_m),Spin(4:6,ipu,i_y,1,ipm),Spin(4:6,i_x,i_y,1,ipm)))*dble(i*j)
!!!!!!!!!!!!!!!!!!
! (100) surface (make a drawing...) the rotation is from v to w
!!!!!!!!!!!!!!!!!!
         charge=charge+(area(Spin(4:6,ipu,i_y,1,i_m),Spin(4:6,ipu,ipv,1,i_m),Spin(4:6,ipu,ipv,1,ipm))+ &
          area(Spin(4:6,ipu,i_y,1,i_m),Spin(4:6,ipu,ipv,1,ipm),Spin(4:6,ipu,i_y,1,ipm)))*dble(i*j)
!!!!!!!!!!!!!!!!!!
! (010) surface (make a drawing...) the rotation is from -u to w
!!!!!!!!!!!!!!!!!!
         charge=charge+(area(Spin(4:6,ipu,ipv,1,i_m),Spin(4:6,i_x,ipv,1,i_m),Spin(4:6,i_x,ipv,1,ipm))+ &
          area(Spin(4:6,ipu,ipv,1,i_m),Spin(4:6,i_x,ipv,1,ipm),Spin(4:6,ipu,ipv,1,ipm)))*dble(i*j)
!!!!!!!!!!!!!!!!!!
! (-100) surface (make a drawing...) the rotation is from w to v
!!!!!!!!!!!!!!!!!!
         charge=charge+(area(Spin(4:6,i_x,i_y,1,i_m),Spin(4:6,i_x,i_y,1,ipm),Spin(4:6,i_x,ipv,1,ipm))+ &
          area(Spin(4:6,i_x,i_y,1,i_m),Spin(4:6,i_x,ipv,1,ipm),Spin(4:6,i_x,ipv,1,i_m)))*dble(i*j)

        enddo
       enddo

       if (charge.gt.1.0d-10) qp=charge/4.0d0/nmag_motif
       if (charge.lt.-1.0d-10) qm=charge/4.0d0/nmag_motif


       end subroutine local_topo_2D_SL

       subroutine local_topo_3D_SL(i_x,i_y,i_z,i_motif,qm,qp,spin,shape_spin,masque,shape_masque)
       use m_parameters, only : Periodic_log
       use m_vector, only : area
       implicit none
! input
       integer, intent(in) :: shape_spin(5),shape_masque(4)
       real(kind=8), intent(in) :: spin(shape_spin(1),shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5))
       integer, intent(in) :: masque(shape_masque(1),shape_masque(2),shape_masque(3),shape_masque(4))
       integer, intent(in) :: i_x,i_y,i_motif,i_z
       real(kind=8), intent(out) :: qm,qp
! internal
       integer :: i,j,ipu,ipv,i_m,ipm
       real(kind=8) :: charge,nmag_motif

       qm=0.0d0
       qp=0.0d0
       charge=0.0d0
       nmag_motif=real(shape_spin(5))

       do j=-1,1,2
        do i=-1,1,2

        if (Periodic_log(1).and.Periodic_log(2)) then
         ipu=mod(i_x+i-1+shape_spin(2),shape_spin(2))+1
         ipv=mod(i_y+j-1+shape_spin(3),shape_spin(3))+1
        else
         if (i_x+i.eq.shape_spin(2)) then
          ipu=i_x
          else
          ipu=i_x+i
         endif
         if (i_y+j.eq.shape_spin(3)) then
          ipv=i_y
          else
          ipv=i_y+1
         endif
        endif

         do i_m=1,shape_spin(5)

          charge=charge+(area(Spin(4:6,i_x,i_y,i_z,i_m),Spin(4:6,ipu,i_y,i_z,i_m),Spin(4:6,ipu,ipv,i_z,i_m))+ &
           area(Spin(4:6,i_x,i_y,i_z,i_m),Spin(4:6,ipu,ipv,i_z,i_m),Spin(4:6,i_x,ipv,i_z,i_m)))*dble(i*j)

         enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! unitl the end, i_m=1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         i_m=1
!!!!!!!!!!!!!!!!!!
! (0-10) surface (make a drawing...) the rotation is from u to w
!!!!!!!!!!!!!!!!!!
         ipm=2

         charge=charge+(area(Spin(4:6,i_x,i_y,i_z,i_m),Spin(4:6,ipu,i_y,i_z,i_m),Spin(4:6,ipu,i_y,i_z,ipm))+ &
          area(Spin(4:6,i_x,i_y,i_z,i_m),Spin(4:6,ipu,i_y,i_z,ipm),Spin(4:6,i_x,i_y,i_z,ipm)))*dble(i*j)
!!!!!!!!!!!!!!!!!!
! (100) surface (make a drawing...) the rotation is from v to w
!!!!!!!!!!!!!!!!!!
         charge=charge+(area(Spin(4:6,ipu,i_y,1,i_m),Spin(4:6,ipu,ipv,1,i_m),Spin(4:6,ipu,ipv,1,ipm))+ &
          area(Spin(4:6,ipu,i_y,i_z,i_m),Spin(4:6,ipu,ipv,i_z,ipm),Spin(4:6,ipu,i_y,i_z,ipm)))*dble(i*j)
!!!!!!!!!!!!!!!!!!
! (010) surface (make a drawing...) the rotation is from -u to w
!!!!!!!!!!!!!!!!!!
         charge=charge+(area(Spin(4:6,ipu,ipv,i_z,i_m),Spin(4:6,i_x,ipv,i_z,i_m),Spin(4:6,i_x,ipv,i_z,ipm))+ &
          area(Spin(4:6,ipu,ipv,i_z,i_m),Spin(4:6,i_x,ipv,i_z,ipm),Spin(4:6,ipu,ipv,i_z,ipm)))*dble(i*j)
!!!!!!!!!!!!!!!!!!
! (-100) surface (make a drawing...) the rotation is from w to v
!!!!!!!!!!!!!!!!!!
         charge=charge+(area(Spin(4:6,i_x,i_y,i_z,i_m),Spin(4:6,i_x,i_y,i_z,ipm),Spin(4:6,i_x,ipv,i_z,ipm))+ &
          area(Spin(4:6,i_x,i_y,i_z,i_m),Spin(4:6,i_x,ipv,i_z,ipm),Spin(4:6,i_x,ipv,i_z,i_m)))*dble(i*j)

        enddo
       enddo

       if (charge.gt.1.0d-10) qp=charge/4.0d0/nmag_motif
       if (charge.lt.-1.0d-10) qm=charge/4.0d0/nmag_motif

       end subroutine local_topo_3D_SL

       subroutine local_vortex_2D(i_x,i_y,v,spin,shape_spin,masque,shape_masque)
       use m_parameters, only : Periodic_log
       use m_vector, only : cross
       implicit none
! input
       integer, intent(in) :: shape_spin(5),shape_masque(4)
       real(kind=8), intent(in) :: spin(shape_spin(1),shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5))
       integer, intent(in) :: masque(shape_masque(1),shape_masque(2),shape_masque(3),shape_masque(4))
       integer, intent(in) :: i_x,i_y
       real(kind=8), intent(out) :: v(3)
! internal
       integer :: i,j,ipu,ipv

       v=0.0d0

       do j=-1,1,2
        do i=-1,1,2

        if (Periodic_log(1).and.Periodic_log(2)) then
         ipu=mod(i_x+i-1+shape_spin(2),shape_spin(2))+1
         ipv=mod(i_y+j-1+shape_spin(3),shape_spin(3))+1
        else
         if (i_x.eq.shape_spin(2)) then
          ipu=i_x
          else
          ipu=i_x+i
         endif
         if (i_y.eq.shape_spin(3)) then
          ipv=i_y
          else
          ipv=i_y+1
         endif
        endif

        v=v+(cross(Spin(4:6,i_x,i_y,1,1),Spin(:,ipu,i_y,1,1))+cross(Spin(4:6,ipu,i_y,1,1),Spin(4:6,ipu,ipv,1,1)) &
        +cross(Spin(4:6,ipu,ipv,1,1),Spin(4:6,i_x,ipv,1,1))+cross(Spin(4:6,i_x,ipv,1,1),Spin(4:6,i_x,i_y,1,1)))*dble(i*j)

        enddo
       enddo

       end subroutine local_vortex_2D

       end module
