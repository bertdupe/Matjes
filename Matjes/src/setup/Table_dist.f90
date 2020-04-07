module m_table_dist
interface Tdist
    module procedure Tdist_basic
    module procedure Tdist_SL
end interface Tdist

interface Tdir
    module procedure Tdir_basic
end interface Tdir

private
public :: get_table_of_distance,Tdir
contains
! Tdist: routine that gives the table of distances up to N nearest neighbors
! Tdir: routine that gives the direction of one of the neighbors up to N nearest neighbors
!




subroutine get_table_of_distance(r,N_Nneigh,world,my_motif,indexNN,n_mag,d)
use m_user_info
use m_derived_types, only : cell
use m_convert
implicit none
real(kind=8), intent(in) :: r(:,:)
type(cell), intent(in) :: my_motif
integer, intent(in) :: N_Nneigh,world(:),indexNN(:,:),n_mag
real(kind=8), intent(inout) :: d(:,:)
! internal variable
integer :: Nei_z,Nei_il,phase,i,shape_table(2),j
real(kind=8) :: time
character(len=50) :: form

Nei_z=0
Nei_il=0
phase=1
time=0.0d0

! calculate distance of NN and table of nearest neighbors
call user_info(6,time,'Calculating the table of distances',.false.)

if ((n_mag.eq.1).and.(phase.eq.1)) then
    call Tdist(r,N_Nneigh,world,my_motif,d(:,1))
else
    call Tdist(r,N_Nneigh,Nei_z,Nei_il,world,phase,my_motif,d(:,:))
endif

call user_info(6,time,'done',.false.)

shape_table=shape(d)
form=convert('(',shape_table(1),'(2x,f10.6))')
write(6,'(a)') ''
write(6,'(a)') 'table of distances'
do i=1,shape_table(2)
   write(6,form) (d(j,i),j=1,shape_table(1))
enddo
write(6,'(a)') ''

end subroutine
!-----------------------------------------
!simple case (and the more commun):
!1 atom per unit cell
!no supercell
!
subroutine Tdist_basic(r,k,N,motif,tabledist)
use m_derived_types, only : cell
use m_vector , only : norm
implicit none
integer, intent(in) :: k,N(:)
!       real (kind=8) :: Tdist(k,phase)
! input variable
real (kind=8), intent(in) :: r(3,3)
type (cell), intent(in) :: motif
!dummy variable
integer :: j,m,pos_min,nmag
real (kind=8) :: u,v,w,minimum,test_vec(3)
! allocation of the output data
real(kind=8), intent(inout) :: tabledist(:)
! dummy results
real(kind=8) :: res((k+1)**size(N)*count(motif%atomic(:)%moment.gt.0.0d0))
!        real(kind=8):: res(36)


nmag=count(motif%atomic(:)%moment.gt.0.0d0)

res=0.0d0
test_vec=0.0d0
j=0
#ifdef CPP_DEBUG
    write(6,*) size(N), "dimension of the system"
    write(6,*) nmag, "magnetic atoms in the motif"
#endif

! 1D case
if (size(N).eq.1) then
   u=0.0d0
   do while (int(u).le.k)
      j=j+1
      res(j)= &
      norm(r(1,:)*(u+motif%atomic(1)%position(1))+r(2,:)*motif%atomic(1)%position(2)+r(3,:)*motif%atomic(1)%position(3))
         u=u+1.0d0
   enddo
!-----------------------------------------
! 2D case
elseif (size(N).eq.2) then
   u=0.0d0
   do while (int(u).le.k)
      v=0.0d0
      do while (int(v).le.k)
         j=j+1
         test_vec=r(1,:)*(u+motif%atomic(1)%position(1))+r(2,:)*(v+motif%atomic(1)%position(2))+r(3,:)*motif%atomic(1)%position(3)
           res(j)=norm(test_vec)
         v=v+1.0d0
      enddo
      u=u+1.0d0
   enddo
else
!-----------------------------------------
! 3D case
   u=0.0d0
   do while (int(u).le.k)
      v=0.0d0
      do while (int(v).le.k)
         w=0.0d0
          do while (int(w).le.k)
             j=j+1
             res(j)= &
          norm(r(1,:)*(u+motif%atomic(1)%position(1))+r(2,:)*(v+motif%atomic(1)%position(2))+r(3,:)*(w+motif%atomic(1)%position(3)))
             w=w+1.0d0
          enddo
         v=v+1.0d0
      enddo
      u=u+1.0d0
   enddo
endif

m=1
minimum=minval(res(:),1,mask=(abs(res).gt.1.0d-8))
tabledist(m)=minimum
if (m.ne.k) then
    do j=1,size(res(:))
       if (abs(res(j)).lt.1.0d-8) cycle
       m=m+1
       pos_min=minloc(res(:),1,mask=(res.gt.(minimum+1.0d-5)))
       if (pos_min.eq.0) exit
       tabledist(m)=res(pos_min)
       minimum=res(pos_min)
       if (m.eq.k) exit
    enddo
endif

#ifdef CPP_DEBUG
      write(6,*) 'table of distance'
      write(6,*) tabledist(:)
#endif

end subroutine Tdist_basic

!-----------------------------------------
!complex case: the superlattice case and the one with more than one atom in the unit cell
!nmag atom per unit cell phase=1
!supperlatice is present: phase=2
!
       subroutine Tdist_SL(r,k,Nei_z,Nei_il,N,phase,motif,tabledist)
       use m_derived_types, only :cell
       use m_vector , only : norm
       implicit none
       integer, intent(in) :: k,N(:),Nei_z,phase,Nei_il
! input variable
       real (kind=8), intent(in) :: r(3,3)
       type (cell), intent(in) :: motif
!dummy variable
       integer :: i,j,l,m,pos_min,nop,nip,nmag,boundary
       real (kind=8) :: u,v,w,minimum
! allocation of the output data
       real(kind=8), intent(inout) :: tabledist(:,:)
! dummy results
       real(kind=8) :: res(max(k+1,Nei_z+1,Nei_il+1)**size(N)*count(motif%i_mom),phase)

       if (phase.ge.2) then
        if (size(N).eq.1) then
         nop=count((sum(motif%pos(:,2:3),2).gt.1.0d-8).and.(motif%i_mom))
         nip=count((motif%pos(:,1).lt.1.0d-8).and.(motif%i_mom))
        else
         nop=count((motif%pos(:,3).gt.1.0d-8).and.(motif%i_mom))
         nip=count((motif%pos(:,3).lt.1.0d-8).and.(motif%i_mom))
        endif
       else
        nmag=count(motif%i_mom)
       endif

       res=0.0d0
       boundary=max(k,Nei_z,Nei_il)

#ifdef CPP_DEBUG
       if (phase.ge.2) then
        write(6,*) "superlattice selected"
        write(6,*) size(N), "dimensional system"
        write(6,*) nop, "magnetic atoms out of world"
        write(6,*) nip, "magnetic atoms in world"
       else
        write(6,*) size(N), "dimension of the system"
        write(6,*) nmag, "magnetic atoms in the motif"
       endif
#endif

! 1D case
       if (size(N).eq.1) then
        u=0.0d0
        do while (int(u).le.boundary)
         if (phase.eq.1) then
          j=0
          do i=1,size(motif%i_mom)
           if (.not.motif%i_mom(i)) cycle
           j=j+1
           res(int(u)*nmag+j,1)= &
            norm(r(1,:)*(u+motif%pos(i,1))+r(2,:)*motif%pos(2,1)+r(2,:)*motif%pos(3,1))
          enddo
          else
           j=0
           l=0
           do i=1,size(motif%i_mom)
            if (.not.motif%i_mom(i)) cycle
            if (abs(sum(motif%pos(i,2:3))).gt.1.0d-8) then
             j=j+1
             res(int(u)*nop+j,2)= &
              norm(r(1,:)*(u+motif%pos(i,1))+r(2,:)*motif%pos(i,2)+r(3,:)*motif%pos(i,3))
            else
             l=l+1
             res(int(u)*nip+l,1)= &
              norm(r(1,:)*(u+motif%pos(i,1))+r(2,:)*motif%pos(i,2)+r(3,:)*motif%pos(i,3))
            endif
           enddo
          endif
         u=u+1.0d0
        enddo
!-----------------------------------------
! 2D case
       elseif (size(N).eq.2) then
        u=0.0d0
        do while (int(u).le.boundary)
        v=0.0d0
         do while (int(v).le.boundary)
          if (phase.eq.1) then
          j=0
           do i=1,size(motif%i_mom)
            if (.not.motif%i_mom(i)) cycle
            j=j+1
            res(int(u)*nmag*max(k,Nei_z)+int(v)*nmag+j,1)= &
            norm(r(1,:)*(u+motif%pos(i,1))+r(2,:)*(v+motif%pos(i,2))+r(3,:)*motif%pos(i,3))
           enddo
           ! phase=2
           else
           j=0
           l=0
           do i=1,size(motif%i_mom)
            if (.not.motif%i_mom(i)) cycle
            if (abs(motif%pos(i,3)).gt.0.0d0) then
             j=j+1
             res(int(u)*nop*boundary+int(v)*nop+j,2)= &
             norm(r(1,:)*(u+motif%pos(i,1))+r(2,:)*(v+motif%pos(i,2))+r(3,:)*motif%pos(i,3))
            else
             l=l+1
             res(int(u)*nip*max(k,Nei_z)+int(v)*nip+l,1)= &
             norm(r(1,:)*(u+motif%pos(i,1))+r(2,:)*(v+motif%pos(i,2))+r(3,:)*motif%pos(i,3))
            endif
           enddo
           endif
          v=v+1.0d0
         enddo
        u=u+1.0d0
        enddo
       else
!-----------------------------------------
! 3D case
        j=0
        l=0
        u=0.0d0
        do while (int(u).le.boundary)
        v=0.0d0
         do while (int(v).le.boundary)
         w=0.0d0
          do while (int(w).le.boundary)
          if (phase.eq.1) then
          j=0
           do i=1,size(motif%i_mom)
            if (.not.motif%i_mom(i)) cycle
            j=j+1
            res(int(u)*nmag*boundary**2+int(v)*nmag*boundary+int(w)*nmag+j,1)= &
            norm(r(1,:)*(u+motif%pos(i,1))+r(2,:)*(v+motif%pos(i,2))+r(3,:)*(w+motif%pos(i,3)))
           enddo
           ! phase=2 or 3
           else
           do i=1,size(motif%i_mom)
            if (.not.motif%i_mom(i)) cycle
            if (abs(motif%pos(i,3)).gt.1.0d-8) then
             j=j+1
             res(j,2)=norm(r(1,:)*(u+motif%pos(i,1))+r(2,:)*(v+motif%pos(i,2))+r(3,:)*(motif%pos(i,3)))
            else
             l=l+1
             res(l,1)=norm(r(1,:)*(u+motif%pos(i,1))+r(2,:)*(v+motif%pos(i,2))+r(3,:)*motif%pos(i,3))

             res(l,3)=norm(r(1,:)*(u+motif%pos(i,1))+r(2,:)*(v+motif%pos(i,2))+r(3,:)*(w+1.0d0+motif%pos(i,3)))
            endif
           enddo
           endif
           w=w+1.0d0
          enddo
          v=v+1.0d0
         enddo
        u=u+1.0d0
        enddo
       endif
       
! routine that finds the minimum

       do i=1,phase
       if (count(res(:,i).eq.0).eq.(k+1)**size(N)*count(motif%i_mom)) cycle
       m=1
       minimum=minval(res(:,i),1,mask=(dabs(res(:,i)).gt.1.0d-8))
       tabledist(m,i)=minimum
       if (m.eq.size(tabledist(:,i))) cycle
        do j=1,size(res(:,i))
         if (abs(res(j,i)).lt.1.0d-8) cycle
         m=m+1
         pos_min=minloc(res(:,i),1,mask=(res(:,i).gt.(1.0d-5+minimum)))
         if (pos_min.eq.0) exit
         tabledist(m,i)=res(pos_min,i)
         minimum=res(pos_min,i)
         if (m.eq.size(tabledist(:,i))) exit
        enddo
       enddo

#ifdef CPP_DEBUG
       do i=1,phase
       write(6,*) 'dimension', i
        do j=1,boundary
        write(6,*) tabledist(j,i)
        enddo
       enddo
#endif

       end subroutine Tdist_SL

!-----------------------------------------
!simple case (and the more commun):
!1 atom per unit cell
!no supercell
!
subroutine Tdir_basic(r,k,N,motif,tabledir)
use m_derived_types, only : cell
use m_vector , only : norm
implicit none
integer, intent(in) :: k,N(:)
!       real (kind=8) :: Tdist(k,phase)
! input variable
real (kind=8), intent(in) :: r(3,3)
type (cell), intent(in) :: motif
!dummy variable
integer :: j,m,pos_min,nmag
real (kind=8) :: u,v,w,minimum,test_vec(3)
! allocation of the output data
real(kind=8), intent(inout) :: tabledir(:,:)
! dummy results
real(kind=8) :: res((k+1)**size(N)*count(motif%i_mom),4)
!        real(kind=8):: res(36)


nmag=count(motif%i_mom)

res=0.0d0
test_vec=0.0d0
j=0
#ifdef CPP_DEBUG
       write(6,*) size(N), "dimension of the system"
       write(6,*) nmag, "magnetic atoms in the motif"
#endif

! 1D case
if (size(N).eq.1) then
   u=0.0d0
   do while (int(u).le.k)
      j=j+1
      res(j,1)= &
      norm(r(1,:)*(u+motif%atomic(1)%position(1))+r(2,:)*motif%atomic(1)%position(2)+r(3,:)*motif%atomic(1)%position(3))
         res(j,2:4)=r(1,:)*(u+motif%atomic(1)%position(1))+r(2,:)*motif%atomic(1)%position(2)+r(3,:)*motif%atomic(1)%position(3)
      u=u+1.0d0
   enddo
!-----------------------------------------
! 2D case
elseif (size(N).eq.2) then
   u=0.0d0
   do while (int(u).le.k)
      v=0.0d0
      do while (int(v).le.k)
         j=j+1
         test_vec=r(1,:)*(u+motif%atomic(1)%position(1))+r(2,:)*(v+motif%atomic(1)%position(2))+r(3,:)*motif%atomic(1)%position(3)
         res(j,1)=norm(test_vec)
         res(j,2:4)=test_vec
         v=v+1.0d0
      enddo
      u=u+1.0d0
   enddo
else
!-----------------------------------------
! 3D case
   u=0.0d0
   do while (int(u).le.k)
      v=0.0d0
      do while (int(v).le.k)
         w=0.0d0
         do while (int(w).le.k)
            j=j+1
            res(j,1)= &
            norm(r(1,:)*(u+motif%atomic(1)%position(1))+r(2,:)*(v+motif%atomic(1)%position(2))+r(3,:)*(w+motif%atomic(1)%position(3)))
              res(j,2:4)=r(1,:)*(u+motif%atomic(1)%position(1))+r(2,:)*(v+motif%atomic(1)%position(2))+r(3,:)*(w+motif%atomic(1)%position(3))
            w=w+1.0d0
         enddo
         v=v+1.0d0
      enddo
      u=u+1.0d0
   enddo
endif

m=1
pos_min=minloc(res(:,1),1,mask=(abs(res(:,1)).gt.1.0d-8))
tabledir(:,m)=res(pos_min,2:4)
minimum=res(pos_min,1)
if (m.ne.k) then
   do j=1,size(res(:,1))
      if (abs(res(j,1)).lt.1.0d-8) cycle
      m=m+1
      pos_min=minloc(res(:,1),1,mask=(res(:,1).gt.(minimum+1.0d-5)))
      if (pos_min.eq.0) exit
      tabledir(:,m)=res(pos_min,2:4)
      minimum=res(pos_min,1)
      if (m.eq.k) exit
   enddo
endif

#ifdef CPP_DEBUG
      write(6,*) 'table of distance'
      write(6,*) tabledir(:,:)
#endif

end subroutine Tdir_basic

end module m_table_dist
