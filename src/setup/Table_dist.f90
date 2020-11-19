module m_table_dist
interface Tdist
    module procedure Tdist_basic
end interface Tdist

private
public :: get_table_of_distance
contains
! Tdist: routine that gives the table of distances up to N nearest neighbors
! Tdir: routine that gives the direction of one of the neighbors up to N nearest neighbors
!




subroutine get_table_of_distance(r,N_Nneigh,world,my_motif,n_mag,d)
use m_user_info
use m_derived_types, only : t_cell
use m_convert
implicit none
real(kind=8), intent(in) :: r(:,:)
type(t_cell), intent(in) :: my_motif
integer, intent(in) :: N_Nneigh,world(:),n_mag
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

call Tdist(r,N_Nneigh,world,my_motif,d(:,1))

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
use m_derived_types, only : t_cell
use m_vector , only : norm
implicit none
integer, intent(in) :: k,N(:)
!       real (kind=8) :: Tdist(k,phase)
! input variable
real (kind=8), intent(in) :: r(3,3)
type(t_cell), intent(in) :: motif
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

end module m_table_dist
