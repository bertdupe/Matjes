module m_indexation
interface numneigh
   module procedure numneigh_simple
end interface numneigh

private
public :: get_num_neighbors
contains

subroutine get_num_neighbors(N_Nneigh,d,r,world,my_motif,indexNN)
use m_derived_types, only : t_cell
implicit none
type(t_cell), intent(in) :: my_motif
integer, intent(in) :: N_Nneigh,world(:)
real(kind=8), intent(in) :: d(:,:),r(:,:)
integer, intent(inout) :: indexNN(:,:)
! internal variables
integer :: phase,Nei_z,Nei_il

phase=1
Nei_z=0
Nei_il=0

call numneigh(N_Nneigh,d(:,1),r,world,my_motif,indexNN(:,1))

end subroutine get_num_neighbors
! calculates the number of neighbours stored in indexNN

subroutine numneigh_simple(k,rad,r,world,motif,indexNN)
use m_derived_types
use m_vector , only : norm
implicit none
! input variable
integer, intent(in) :: k,world(:)
real (kind=8), intent(in) :: rad(:),r(3,3)
type(t_cell), intent(in) :: motif
! variable that goes out
integer, intent(inout) :: indexNN(:)
!dummy variable
integer :: i,j,size_mag
real (kind=8) :: u,v,w,val

indexNN=0
size_mag=count(motif%atomic(:)%moment.gt.0.0d0)

if (size(world).eq.1) then
!----------------------------------------
! cas 1D
    do i=1,k
       u=-dble(i)
         do while (int(dabs(u)).le.i)
           do j=1,size_mag
             if (motif%atomic(j)%moment.lt.1.0d-8) cycle
             val=norm(r(1,:)*(u+motif%atomic(j)%position(1))+r(2,:)*motif%atomic(j)%position(2)+ &
         & r(3,:)*motif%atomic(j)%position(3))
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
                 if (motif%atomic(j)%moment.lt.1.0d-8) cycle
                 val=norm(r(1,:)*(u+motif%atomic(j)%position(1))+r(2,:)*(v+motif%atomic(j)%position(2))+ &
          & r(3,:)*motif%atomic(j)%position(3))
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
                      if (motif%atomic(j)%moment.lt.1.0d-8) cycle
                      val=norm(r(1,:)*(u+motif%atomic(j)%position(1))+r(2,:)*(v+motif%atomic(j)%position(2))+ &
          & r(3,:)*(w+motif%atomic(j)%position(3)))
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

end module m_indexation
