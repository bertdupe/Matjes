module m_vector
implicit none


   interface project
       module procedure vec_project
   end interface 

   interface normalize
       module procedure vec_normalize
       module procedure vec_normalize_eps
   end interface

   interface norm
       module procedure norm_real
       module procedure norm_int
   end interface norm

   interface cross
       module procedure cross_real_simple
       module procedure cross_real
!       module procedure cross_realp
       module procedure cross_int
   end interface cross

   interface norm_cross
       module procedure norm_cross_real
       module procedure norm_cross_int
   end interface norm_cross

   interface distance
       module procedure d_1d
       module procedure d_2d
   end interface distance

   interface min_dist
       module procedure one_d_dist
       module procedure min_dist_2D
       module procedure min_dist_3D
   end interface min_dist

   interface calc_ang
      module procedure calc_ang_real
      module procedure calc_ang_realarr
   end interface calc_ang

   interface calculate_damping
      module procedure damping_2V
   end interface

public

contains

! ==============================================================
!> projects vector vec on reference
pure subroutine vec_project(vec,ref)
    real(8),intent(inout)   :: vec(:,:)
    real(8),intent(in)      :: ref(:,:)
    !internal
    real(8)     :: prod
    integer     :: i

    do i=1,size(vec,2)
        prod=DOT_PRODUCT(vec(:,i),ref(:,i))
        vec(:,i)=vec(:,i)-prod*ref(:,i)
    enddo
end subroutine

! ==============================================================
!> Normalizes a vector
pure subroutine vec_normalize(vec)
    real(8),intent(inout)   :: vec(:,:)
    !internal
    real(8)     :: nrm(size(vec,2))
    integer     :: i

    nrm=norm2(vec,1)
    forall(i=1:size(nrm), nrm(i)>1.0d-8 ) vec(:,i)=vec(:,i)/nrm(i)
end subroutine


pure subroutine vec_normalize_eps(vec,eps)
    real(8),intent(inout)   :: vec(:,:)
    real(8),intent(in)      :: eps
    !internal
    real(8)     :: nrm(size(vec,2))
    integer     :: i

    nrm=norm2(vec,1)
    forall(i=1:size(nrm), nrm(i)>eps ) vec(:,i)=vec(:,i)/nrm(i)
end subroutine



! ==============================================================
!> Calculate angle between two 3-vectors
real(kind=8) function calc_ang_real(n1,n2)
implicit none
real(kind=8), intent(in) :: n1(:), n2(:) !n1 and n2 have to be normalized
!real(kind=8) :: calc_ang
!internal variables
real(kind=8) :: prod,tmp
integer :: N

N=size(n1)

tmp = norm_cross(n1,n2,1,N)
prod = dot_product(n1,n2)

calc_ang_real = datan2(tmp,prod)
end function calc_ang_real


function calc_ang_realarr(n1,n2)result(res)
real(8), intent(in),contiguous,target :: n1(:,:), n2(:,:) !n1 and n2 have to be normalized
!real(kind=8) :: calc_ang
!internal variables
real(8),pointer     :: flat_1(:),flat_2(:),flat_tmp(:)
real(8)             ::  res(size(n1,2))
real(8),target      :: tmp(3,size(n1,2))
real(8)             :: nrm(size(n1,2))
real(8)             :: prod(size(n1,2))
integer :: N

N=size(n1)
flat_1(1:N)=>n1
flat_2(1:N)=>n2
flat_tmp(1:N)=>tmp

flat_tmp = cross(flat_1,flat_2,1,N)
nrm=norm2(tmp,1)
tmp=n1*n2
prod=sum(tmp,1)

res = datan2(nrm,prod)
end function 


! ==============================================================

function damping_2V(A,B)
! calculates Ax(AxB)=-(A.(A.B)-B)
implicit none
real(kind=8), intent(in) :: A(:),B(:)
real(kind=8), dimension(size(A)) :: damping_2V
! internal
integer :: i
real(kind=8) :: dummy

damping_2V=0.0d0
do i=1,size(A)/3
  dummy=dot_product( A((i-1)*3+1:i*3) ,B((i-1)*3+1:i*3) )
  damping_2V((i-1)*3+1:i*3) = -(A((i-1)*3+1:i*3) * dummy - B((i-1)*3+1:i*3))
enddo

end function damping_2V

! ==============================================================

        real(kind=8) function one_d_dist(nn,ii,jj)
        ! calculates the distance between i and j in one Dimension
        implicit none
        real(kind=8), intent(in) :: nn,ii,jj
        real(kind=8) :: aa,bb,diss
       !do j=1,N
        aa=abs(jj-ii)
        if(ii<jj) then
                bb=abs(nn-jj+ii)
         else
                bb=abs(nn-ii+jj)
        end if
        if (aa<bb) then
                diss = aa
        else
                diss = bb
        end if
        !end do
        one_d_dist=diss
        return
        end function one_d_dist

! ==============================================================

        function min_dist_2D(nx,ny,S1,S2,rx,ry)
        ! calculates the distance between i and j
        implicit none
        real(kind=8), intent(in) :: nx,ny,S1(3),S2(3),rx(3),ry(3)
        real(kind=8),dimension(3) :: min_dist_2D
! internals
        real(kind=8) :: test(5),vector(5,3)
        integer :: pos,j

         min_dist_2D(:)=0.0d0

        vector(1,:)=S2-S1
        vector(2,:)=S2-S1-nx*rx
        vector(3,:)=S2-S1-nx*rx-ny*ry
        vector(4,:)=S2-S1-ny*ry
        vector(5,:)=S2-S1-nx*rx+ny*ry

        do j=1,5
         test(j)=sqrt(vector(j,1)**2+vector(j,2)**2+vector(j,3)**2)
        enddo

        pos=minloc(test,1)
        min_dist_2D=vector(pos,:)
        end function min_dist_2D

! ==============================================================

        function min_dist_3D(nx,ny,nz,S1,S2,rx,ry,rz)
        ! calculates the distance between i and j
        implicit none
        real(kind=8) :: min_dist_3D(3)
        real(kind=8), intent(in) :: nx,ny,nz,S1(3),S2(3),rx(3),ry(3),rz(3)
! internals
        real(kind=8) :: test(11),vector(11,3)
        integer :: pos,j

        vector(1,:)=S2-S1
        vector(2,:)=S2-S1-nx*rx
        vector(3,:)=S2-S1-nx*rx-ny*ry
        vector(4,:)=S2-S1-ny*ry
        vector(5,:)=S2-S1-nx*rx+ny*ry

        vector(6,:)=S2-S1-nx*rx-nz*rz
        vector(7,:)=S2-S1-nx*rx-ny*ry-nz*rz
        vector(8,:)=S2-S1-ny*ry-nz*rz
        vector(9,:)=S2-S1-nx*rx+ny*ry-nz*rz

        vector(10,:)=S2-S1-nx*rx+nz*rz
        vector(11,:)=S2-S1-ny*ry+nz*rz

        do j=1,11
         test(j)=sqrt(vector(j,1)**2+vector(j,2)**2+vector(j,3)**2)
        enddo

        pos=minloc(test,1)
        min_dist_3D=vector(pos,:)
        end function min_dist_3D

! calculate the phi angle
      real(kind=8) function phi(a,x,y,dim_lat,net)
      integer, intent(in) :: dim_lat(3)
      real(kind=8), intent(in) :: net(3,3),a(2)
      real(kind=8), intent(in) :: x,y
!internal
      real(kind=8) :: u,v,d,test,test_sin

      d=sqrt((a(1)-x)**2+(a(2)-y)**2)
      if (abs(d).gt.1.0d-8) then
       phi=acos((a(1)-x)/d)
       test_sin=asin((a(2)-y)/d)
       if (test_sin.lt.0.0d0) phi=-phi+2.0d0*3.14159265359d0
       else
       phi=0.0d0
       return
      endif

      u=-1.0d0
      do while (u.le.1.0d0)
       v=-1.0d0
       do while (v.le.1.0d0)
        test=sqrt((a(1)+u*real(dim_lat(1),8)*net(1,1)+v*real(dim_lat(2),8)*net(2,1)-x)**2+ &
     & (a(2)+u*real(dim_lat(1),8)*net(1,2)+v*real(dim_lat(2),8)*net(2,2)-y)**2)
        if (test.lt.d) then
         d=test
         phi=acos((a(1)+u*real(dim_lat(1),8)*net(1,1)+v*real(dim_lat(2),8)*net(2,1)-x)/test)
         test_sin=asin((a(2)+u*real(dim_lat(1),8)*net(1,2)+v*real(dim_lat(2),8)*net(2,2)-y)/test)
         if (test_sin.lt.0.0d0) phi=-phi+2.0d0*3.14159265359d0
        endif
        v=v+1.0d0
       enddo
       u=u+1.0d0
      enddo

      end function phi

! function that gives the shortest distance between two points in the lattice
      real(kind=8) function d_1d(a,x,dim_lat,net)
      integer, intent(in) :: dim_lat(:)
      real(kind=8), intent(in) :: net(:,:),a(:)
      real(kind=8), intent(in) :: x
!internal
      d_1d=minval((/abs(a(1)-x),abs(a(1)+dble(dim_lat(1))*net(1,1)-x),abs(a(1)-dble(dim_lat(1))*net(1,1)-x)/))
      end function d_1d

      real(kind=8) function d_2d(a,x,y,dim_lat,net)
      integer, intent(in) :: dim_lat(:)
      real(kind=8), intent(in) :: net(:,:),a(:)
      real(kind=8), intent(in) :: x,y
!internal
      real(kind=8) :: u,v,test

      d_2d=sqrt((a(1)-x)**2+(a(2)-y)**2)

      u=-1.0d0
      do while (u.le.1.0d0)
       v=-1.0d0
       do while (v.le.1.0d0)
        test=sqrt((a(1)+u*dble(dim_lat(1))*net(1,1)+v*dble(dim_lat(2))*net(2,1)-x)**2+ &
     & (a(2)+u*dble(dim_lat(1))*net(1,2)+v*dble(dim_lat(2))*net(2,2)-y)**2)
        if (test.lt.d_2d) d_2d=test
        v=v+1.0d0
       enddo
       u=u+1.0d0
      enddo

      end function d_2d

! function that inverts a linear 3D system A(3)*X(3,3)=B(3). the lattice vector in X have to be given in column.
! at the end A=AXeqB
      function AXeqB(X,B)
      implicit none
      integer :: AXeqB(3)
      real(kind=8), intent(in) :: X(3,3),B(3)
! dummies
      real(kind=8) :: denom

      denom=M33DET(X)
      if (abs(denom).lt.1.0d-8) then
       write(*,*) "problem in AXeqB"
       stop
      endif

      AXeqB(1)=nint(-(-B(3)*X(2,2)*X(1,3)+B(2)*X(3,2)*X(1,3)+B(3)*X(1,2)*X(2,3)-B(1)*X(3,2)*X(2,3)-B(2)*X(1,2)*X(3,3) &
     & +B(1)*X(2,2)*X(3,3))/denom)

      AXeqB(2)=nint(-(B(3)*X(2,1)*X(1,3)-B(2)*X(3,1)*X(1,3)-B(3)*X(1,1)*X(2,3)+B(1)*X(3,1)*X(2,3)+B(2)*X(1,1)*X(3,3) &
     & -B(1)*X(2,1)*X(3,3))/denom)

      AXeqB(3)=nint((B(3)*X(2,1)*X(1,2)-B(2)*X(3,1)*X(1,2)-B(3)*X(1,1)*X(2,2)+B(1)*X(3,1)*X(2,2)+B(2)*X(1,1)*X(3,2) &
     & -B(1)*X(2,1)*X(3,2))/denom)

      end function









!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real(kind=8) function norm_real(v)
implicit none
real(kind=8), intent(in) :: v(:)
norm_real=sqrt(sum(v**2))
end function norm_real

real(kind=8) function norm_int(v)
implicit none
integer, intent(in) :: v(:)
norm_int=sqrt(sum(dble(v)**2))
end function norm_int







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! work only for vectors of dimension 3N
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real(kind=8) FUNCTION norm_cross_real(a, b,P,N)
implicit none
integer, intent(in) :: N,P
real(kind=8), INTENT(IN) :: a(:), b(:)
real(kind=8) :: cross_real(P:N)
! internal
integer :: i,j

if (mod(N-P+1,3).ne.0) STOP 'error in norm_cross_real'
norm_cross_real=0.0d0

do i=1,(N-P+1)/3
   j=3*(i-1)+P
   cross_real(j)   =  a(j+1) * b(j+2) - a(j+2) * b(j+1)
   cross_real(j+1) =  a(j+2) * b(j)   - a(j)   * b(j+2)
   cross_real(j+2) =  a(j)   * b(j+1) - a(j+1) * b(j)

   norm_cross_real=norm_cross_real+norm(cross_real(j:j+2))**2
enddo

norm_cross_real=sqrt(norm_cross_real)

END FUNCTION norm_cross_real

real(kind=8) FUNCTION norm_cross_int(a, b,P,N)
implicit none
integer, intent(in) :: N,P
integer, INTENT(IN) :: a(:), b(:)
real(kind=8) :: cross_real(P:N)
! internal
integer :: i,j

if (mod(N-P+1,3).ne.0) STOP 'error in norm_cross_int'

norm_cross_int=0.0d0
do i=1,(N-P+1)/3
   j=3*(i-1)+P
   cross_real(j)   = dble( a(j+1) * b(j+2) - a(j+2) * b(j+1) )
   cross_real(j+1) = dble( a(j+2) * b(j)   - a(j)   * b(j+2) )
   cross_real(j+2) = dble( a(j)   * b(j+1) - a(j+1) * b(j)   )

   norm_cross_int=norm_cross_int+norm(cross_real(j:j+2))**2
enddo

norm_cross_int=sqrt(norm_cross_int)

END FUNCTION norm_cross_int



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION cross_real_simple(a, b)
implicit none
real(kind=8), INTENT(IN) :: a(3), b(3)
real(kind=8) :: cross_real_simple(3)
! internal

cross_real_simple(1) = a(2) * b(3) - a(3) * b(2)
cross_real_simple(2) = a(3) * b(1) - a(1) * b(3)
cross_real_simple(3) = a(1) * b(2) - a(2) * b(1)

END FUNCTION cross_real_simple

!FUNCTION cross_realp(a, b, P, N)result(res)
!    integer, intent(in) :: N,P
!    real(8),pointer, INTENT(IN) :: a(:), b(:)
!    real(8) :: res(P:N)
!
!    res=cross_real(a,b,P,N)
!end function

FUNCTION cross_real(a, b, P, N)
implicit none
integer, intent(in) :: N,P
real(kind=8), INTENT(IN) :: a(:), b(:)
real(kind=8) :: cross_real(P:N)
! internal
integer :: i,j

if (mod(N-P+1,3).ne.0) STOP 'error in cross_real'

do i=1,(N-P+1)/3
   j=3*(i-1)+P
   cross_real(j)   = a(j+1) * b(j+2) - a(j+2) * b(j+1)
   cross_real(j+1) = a(j+2) * b(j)   - a(j)   * b(j+2)
   cross_real(j+2) = a(j)   * b(j+1) - a(j+1) * b(j)
enddo

END FUNCTION cross_real

FUNCTION cross_int(a, b, P, N)
implicit none
integer, intent(in) :: N,P
integer, INTENT(IN) :: a(:), b(:)
real(kind=8) :: cross_int(P:N)
! internal
integer :: i,j

if (mod(N-P+1,3).ne.0) STOP 'error in cross_real'

do i=1,(N-P+1)/3
   j=3*(i-1)+P
   cross_int(j)   = dble( a(j+1) * b(j+2) - a(j+2) * b(j+1) )
   cross_int(j+1) = dble( a(j+2) * b(j)   - a(j)   * b(j+2) )
   cross_int(j+2) = dble( a(j)   * b(j+1) - a(j+1) * b(j)   )
enddo

END FUNCTION cross_int



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real(kind=8) FUNCTION M33DET (A)
IMPLICIT NONE
real(kind=8), DIMENSION(3,3), INTENT(IN)  :: A
real(kind=8) :: DET

DET =   A(1,1)*A(2,2)*A(3,3) &
       - A(1,1)*A(2,3)*A(3,2) &
       - A(1,2)*A(2,1)*A(3,3) &
       + A(1,2)*A(2,3)*A(3,1) &
       + A(1,3)*A(2,1)*A(3,2) &
       - A(1,3)*A(2,2)*A(3,1)

M33DET=DET
END FUNCTION M33DET

! function used to calculate the topological charge
real(kind=8) function area(u0,v0,w0)
IMPLICIT NONE
real(kind=8), DIMENSION(3), intent(in) :: u0,v0,w0
real(kind=8), DIMENSION(3,3) :: JACNN
real(kind=8)  :: pr,n_u0,n_v0,n_w0

area=0.0d0
JACNN=0.0d0

n_u0=sqrt(u0(1)**2+u0(2)**2+u0(3)**2)
n_v0=sqrt(v0(1)**2+v0(2)**2+v0(3)**2)
n_w0=sqrt(w0(1)**2+w0(2)**2+w0(3)**2)

if (n_u0.gt.1.0d-8) JACNN(:,1)=u0/n_u0
if (n_v0.gt.1.0d-8) JACNN(:,2)=v0/n_v0
if (n_w0.gt.1.0d-8) JACNN(:,3)=w0/n_w0

pr=1.0d0+dot_product(v0,w0)+dot_product(u0,w0)+dot_product(u0,v0)
if (abs(pr).lt.1.0d-8) then
    area=acos(-1.0d0)*sign(1.0d0,pr)*sign(1.0d0,M33DET(JACNN))
else
    area=2.0d0*atan2(M33DET(JACNN),pr)
endif

end function area

! function used to TEST the calculation of the topological charge
!      real(kind=8) function area_test(u0,v0,w0)
!      IMPLICIT NONE
!      real(kind=8), DIMENSION(3), intent(in) :: u0,v0,w0
!      real(kind=8), DIMENSION(3,3) :: JACNN
!      real(kind=8)  :: pr
!       JACNN(1,1:3)=u0(1:3)
!       JACNN(2,1:3)=v0(1:3)
!       JACNN(3,1:3)=w0(1:3)
!       pr=1.0d0+dot_product(v0,w0)+dot_product(u0,w0)+dot_product(u0,v0)
!       if (abs(pr).lt.1.0d-8) then
!        area_test=acos(-1.0d0)*sign(1.0d0,pr)*sign(1.0d0,M33DET(JACNN))
!       else
!        area_test=2.0d0*atan2(M33DET(JACNN),pr)
!       endif
!      end function area_test

function sorien(u0,v0,w0)
IMPLICIT NONE
real(kind=8) :: sorien(3)
real(kind=8), DIMENSION(3), intent(in) :: u0,v0,w0
real(kind=8), DIMENSION(3,3) :: JACNN
real(kind=8) :: pr
real(kind=8) :: dumy(3)

dumy=u0+v0+w0
JACNN(1,1:3)=u0(:)
JACNN(2,1:3)=v0(:)
JACNN(3,1:3)=w0(:)
pr=1.0d0+dot_product(v0,w0)+dot_product(u0,w0)+dot_product(u0,v0)
if (abs(pr).lt.1.0d-8) pr=1.0d0
sorien=2.0d0*atan2(M33DET(JACNN),pr)*(u0+v0+w0)/sqrt(dumy(1)**2+dumy(2)**2+dumy(3)**2)
end function sorien

! ===============================================================
real(kind=8) function TripleProduct(U,S,T)
Implicit none
!     Vectors, and result Vector
real(kind=8), Dimension(1:3) :: S, T, U

!      res=0.0d0
!      res(1)=(S(2)*T(3)-S(3)*T(2))*U(1)
!      res(2)=(S(3)*T(1)-S(1)*T(3))*U(2)
!      res(3)=(S(1)*T(2)-S(2)*T(1))*U(3)

TripleProduct=(S(2)*T(3)-S(3)*T(2))*U(1)+(S(3)*T(1)-S(1)*T(3))*U(2)+(S(1)*T(2)-S(2)*T(1))*U(3)

end function TripleProduct

end module m_vector
