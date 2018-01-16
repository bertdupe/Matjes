       module m_indexation
       interface numneigh
        module procedure numneigh_simple
        module procedure numneigh_SL
       end interface numneigh
       contains

! calculates the number of neighbours stored in indexNN

       subroutine numneigh_simple(k,rad,r,world,motif,indexNN)
       use m_derived_types
       use m_vector , only : norm
       implicit none
! input variable
       integer, intent(in) :: k,world(:)
       real (kind=8), intent(in) :: rad(:),r(3,3)
       type (cell), intent(in) :: motif
! variable that goes out
       integer, intent(inout) :: indexNN(:)
!dummy variable
       integer :: i,j,size_mag
       real (kind=8) :: u,v,w,val

       indexNN=0
       size_mag=size(motif%i_mom)

       if (size(world).eq.1) then
!----------------------------------------
! cas 1D
       do i=1,k
        u=-dble(i)
        do while (int(dabs(u)).le.i)
         do j=1,size_mag
          if (.not.motif%i_mom(j)) cycle
          val=norm(r(1,:)*(u+motif%pos(j,1))+r(2,:)*motif%pos(j,2)+r(3,:)*motif%pos(j,3))
          if (abs(rad(i)-val).lt.1.0d-8) indexNN(i)=indexNN(i)+1
         enddo
         u=u+1.0d0
        enddo
       enddo
!----------------------------------------
! cas 2D
       elseif (size(world).eq.2) then
       do i=1,k
        u=-dble(i)
        do while (int(dabs(u)).le.i)
         v=-dble(i)
         do while (int(dabs(v)).le.i)
          do j=1,size_mag
           if (.not.motif%i_mom(j)) cycle
           val=norm(r(1,:)*(u+motif%pos(j,1))+r(2,:)*(v+motif%pos(j,2))+r(3,:)*motif%pos(j,3))
           if (dabs(rad(i)-val).lt.1.0d-8) indexNN(i)=indexNN(i)+1
          enddo
          v=v+1.0d0
         enddo
         u=u+1.0d0
        enddo
       enddo
!----------------------------------------
! cas 3D
       else
       do i=1,k
        u=-dble(i)
        do while (int(dabs(u)).le.i)
         v=-dble(i)
         do while (int(dabs(v)).le.i)
          w=-dble(i)
          do while (int(dabs(w)).le.i)
           do j=1,size_mag
            if (.not.motif%i_mom(j)) cycle
            val=norm(r(1,:)*(u+motif%pos(j,2))+r(2,:)*(v+motif%pos(j,2))+r(3,:)*(w+motif%pos(j,3)))
            if (dabs(rad(i)-val).lt.1.0d-8) indexNN(i)=indexNN(i)+1
           enddo
           w=w+1.0d0
          enddo
          v=v+1.0d0
         enddo
         u=u+1.0d0
        enddo
       enddo
       endif

#ifdef CPP_DEBUG
       write(6,*) "simple case selected: no SL"
       write(6,*) "indexNN  :", indexNN(:)
#endif
       end subroutine numneigh_simple

! calculates the number of neighbours stored in indexNN

       subroutine numneigh_SL(k,rad,r,Nei_z,Nei_il,world,motif,phase,indexNN)
       use m_derived_types
       use m_vector , only : norm
       implicit none
! input variable
       integer, intent(in) :: k,Nei_z,world(:),phase,Nei_il
       real (kind=8), intent(in) :: rad(max(k,Nei_z,Nei_il),phase),r(3,3)
       type (cell), intent(in) :: motif
! variable that goes out
       integer, intent(inout) :: indexNN(max(k,Nei_z),phase)
!dummy variable
       integer :: i,j,size_mag
       real (kind=8) :: u,v,w,val

       indexNN=0
       size_mag=size(motif%i_mom)

       if (size(world).eq.1) then
!----------------------------------------
! cas 1D
       do i=1,max(k,Nei_il)
        u=-dble(i)
        do while (int(dabs(u)).le.i)
         do j=1,size_mag
          if (.not.motif%i_mom(j)) cycle
           val=norm(r(1,:)*(u+motif%pos(j,1))+r(2,:)*motif%pos(j,2)+r(3,:)*motif%pos(j,3))
           if ((abs(rad(i,1)-val).lt.1.0d-8).and.(abs(sum(motif%pos(j,2:3))).lt.1.0d-8)) then
            indexNN(i,1)=indexNN(i,1)+1
           elseif (abs(rad(i,2)-val).lt.1.0d-8) then
            indexNN(i,2)=indexNN(i,2)+1
           endif
         enddo
         u=u+1.0d0
        enddo
       enddo
!----------------------------------------
! cas 2D
       elseif (size(world).eq.2) then
       do i=1,max(k,Nei_il)
        u=-dble(i)
        do while (int(dabs(u)).le.i)
         v=-dble(i)
         do while (int(dabs(v)).le.i)
          do j=1,size_mag
           if (.not.motif%i_mom(j)) cycle
            val=norm(r(1,:)*(u+motif%pos(j,1))+r(2,:)*(v+motif%pos(j,2))+r(3,:)*motif%pos(j,3))
            if (i.le.k) then
             if ((dabs(rad(i,1)-val).lt.1.0d-8).and.(abs(motif%pos(j,3)).lt.1.0d-8)) indexNN(i,1)=indexNN(i,1)+1
            endif
            if (i.le.Nei_il) then
             if ((dabs(rad(i,2)-val).lt.1.0d-8).and.(abs(motif%pos(j,3)).gt.1.0d-8)) indexNN(i,2)=indexNN(i,2)+1
            endif
          enddo
          v=v+1.0d0
         enddo
         u=u+1.0d0
        enddo
       enddo
!----------------------------------------
! cas 3D
       elseif (size(world).eq.3) then
       do i=1,max(k,Nei_z,Nei_il)
        u=-dble(i)
        do while (int(dabs(u)).le.i)
         v=-dble(i)
         do while (int(dabs(v)).le.i)
          w=-dble(i)
          do while (int(dabs(w)).le.i)
           do j=1,size_mag
            if (.not.motif%i_mom(j)) cycle
            val=norm(r(1,:)*(u+motif%pos(j,1))+r(2,:)*(v+motif%pos(j,2))+r(3,:)*(w+motif%pos(j,3)))

            if ((abs(rad(i,1)-val).lt.1.0d-8).and.(abs(motif%pos(j,3)).lt.1.0d-8).and.(i.le.k)) indexNN(i,1)=indexNN(i,1)+1

            if ((abs(rad(i,2)-val).lt.1.0d-8).and.(abs(motif%pos(j,3)).gt.1.0d-8).and.(i.le.Nei_il)) indexNN(i,2)=indexNN(i,2)+1

            if ((abs(rad(i,3)-val).lt.1.0d-8).and.(i.le.Nei_z)) indexNN(i,3)=indexNN(i,3)+1
           enddo
           w=w+1.0d0
          enddo
          v=v+1.0d0
         enddo
         u=u+1.0d0
        enddo
       enddo
       endif

#ifdef CPP_DEBUG
       write(6,*) "SuperLattice case selected"
       do i=1,phase
        write(6,*) 'indexNN in dimension  ', i,":", indexNN(:,i)
       enddo
#endif
       end subroutine numneigh_SL

       end module m_indexation
