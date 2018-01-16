      module m_sym_utils
      integer, allocatable :: corners(:,:)
      interface rot_mat
       module procedure rot_mat_I
       module procedure rot_mat_R
      end interface rot_mat
      interface pos_nei
       module procedure pos_nei_2D
       module procedure pos_nei_3D_SL
       module procedure pos_nei_2D_real
       module procedure pos_nei_2D_user
       module procedure pos_between2at_2D
      end interface pos_nei
      contains
!calculate the order of the z axis
      integer function order_zaxis(r)
      use m_constants
      use m_vector, only : norm
      implicit none
      real(kind=8), intent(in) :: r(3,3)
! dummy
      real(kind=8) :: ss

      ss=dot_product(r(1,:),r(2,:))/norm(r(1,:))/norm(r(2,:))

      order_zaxis=360/nint(180.0d0/pi(1.0d0)*dacos(ss))

      order_zaxis=order_zaxis*int(dsign(1.0d0,datan(r(2,2)/r(1,1))))

      end function order_zaxis

! find oriented angles between 2 vectors to go from u to v
      real(kind=8) function angle_oriented(u,v)
      implicit none
      real(kind=8), dimension(3), intent(in) :: u,v
! dumy
      real(kind=8) :: angle,ss

      ss=(u(1)*v(1)+u(2)*v(2)+u(3)*v(3))/(dsqrt(u(1)**2+u(2)**2+u(3)**2)* &
       dsqrt(v(1)**2+v(2)**2+v(3)**2))

      angle=dacos(ss)
      angle_oriented=angle*dsign(1.0d0,dasin(dsqrt(1-ss**2)))

      end function angle_oriented

! ===============================================================
      ! calculate the reciprocal-mesh
      ! written by Lukas Deuchler
      ! date 11/19/2013
      ! email deuchlerl@gmail.com
      !
      integer function mesh_bz(order_z_axis,r,N)
      use m_constants
      use m_vector, only : norm
! intent(in)
      integer, intent(in) :: order_z_axis,N
      real(kind=8), intent(in) :: r(3,3)
!value of the function
      real(kind=8) :: b1(3),b2(3),v1(3),v2(3)
      real(kind=8) :: rot_ang,rotation(3,3),s
      integer :: i,j
      real(kind=8) :: len_v1

      mesh_bz=0

      b1=r(1,:)/2.0d0
      b2=r(2,:)/2.0d0


      OPEN(1,FILE='kpoints',action='write',status='unknown', &
       form='formatted')

      write(1,'(3(f10.7,3x))') 0.0d0,0.0d0,0.0d0

      if (mod(order_z_axis,3).eq.0) then
       rot_ang=30.0d0/360.0d0*pi(2.0d0) 
       rotation=reshape((/dcos(rot_ang),-dsin(rot_ang),0.0d0,&
       dsin(rot_ang),dcos(rot_ang),0.0d0,&
       0.0d0,0.0d0,1.0d0/),(/3,3/))
       v1 = matmul(b1,rotation)
       len_v1 = norm(b1)/dcos(rot_ang)
       s=norm(v1)
       v1 = len_v1*v1/s

       rot_ang=-60.0d0/360.0d0*pi(2.0d0)
       rotation=reshape((/dcos(rot_ang),-dsin(rot_ang),0.0d0,&
       dsin(rot_ang),dcos(rot_ang),0.0d0,&
       0.0d0,0.0d0,1.0d0/),(/3,3/))
       v2 = matmul(v1,rotation)
       s=norm(v2)
       v2 = len_v1*v2/s


       if (dabs(dot_product(b1,v1-v2)).lt.1.0d-7) mesh_bz=1
!per line

       do i=1,N+1
        do j=2,i
         write(1,'(3(f10.7,3x))') v2*dble(i-1)/dble(N)+(v1-v2)*dble(j-1)/dble(N) 
        end do 
       end do      

      else 

       do i=1,N
        do j=2,N
         write(1,'(3(f10.7,3x))') b1*dble(i-1)/dble(N)+b2*dble(j-1)/dble(N)
        end do
       end do

      end if

      close(1)

      if (mesh_bz.eq.0) then
       write(*,'(a)') 'problem building the mesh'
       stop
      endif

      end function mesh_bz

!
! function that gives the position of closest symmetry equivalent non magnetic atom
! between 2 magnetic atoms
!

      function pos_between2at_2D(pos_ref,pos,pos2,r)
      use m_vector, only : norm
      implicit none
      real(kind=8) :: pos_between2at_2D(3)
      real(kind=8), intent(in) :: pos_ref(:),pos(:),pos2(:)
      real(kind=8), intent(in) :: r(:,:)
      !dummy
      integer :: i,j,k
      real(kind=8) :: table_dist(10,2),table_pos(10,3),vec_1(3),vec_2(3)

      pos_between2at_2D=0.0d0

      k=1
      table_dist(1,1)=norm(pos_ref-pos)
      table_dist(1,2)=norm(pos_ref-pos2)
      table_pos(1,:)=pos_ref-pos
      k=2

      do i=-1,1,1
       do j=-1,1,1
        vec_1=pos_ref-pos+dble(i)*r(1,:)+dble(j)*r(2,:)
        vec_2=pos_ref-pos2+dble(i)*r(1,:)+dble(j)*r(2,:)
        table_dist(k,1)=norm(vec_1)
        table_dist(k,2)=norm(vec_2)
        table_pos(k,:)=pos_ref+dble(i)*r(1,:)+dble(j)*r(2,:)
        k=k+1
       enddo
      enddo

      do k=1,10
       if (abs(table_dist(k,1)-table_dist(k,2)).lt.1.0d-4) then
        pos_between2at_2D=table_pos(k,:)
        return
       endif
       if (k.eq.10) then
        write(6,'(a)') ' '
        write(6,'(a)') ' '
        write(6,'(a)') '!!!!! WARNING !!!!!'
        write(6,'(a)') 'No Ir was find at equidistance between the 2 magnetic neighbors'
        write(6,'(a)') 'DM will be purely in plane then'
        write(6,'(a)') ' '
        write(6,'(a)') ' '
        pos_between2at_2D=pos2/2.0d0
        pos_between2at_2D(3)=-1000.0d0
       endif
      enddo

      end function pos_between2at_2D

!
! function that gives the position closest neighbor among the different supercell image
! with periodic boundary condition
!
      function pos_nei_2D(ix,iy,pos,r,dim_lat)
      use m_vector, only : norm
      implicit none
      real(kind=8) :: pos_nei_2D(3)
      integer, intent(in) :: dim_lat(3),pos(3),ix,iy
      real(kind=8), intent(in) :: r(3,3)
      !dummy
      integer :: i,j,k,min_dist
      real(kind=8) :: table_dist(10),table_pos(10,3),ref(3)

      ref=dble(ix-1)*r(1,:)+dble(iy-1)*r(2,:)+0.0d0*r(3,:)
      k=1
      table_dist(1)=norm(dble(pos(1)-1)*r(1,:)+dble(pos(2)-1)*r(2,:)+dble(pos(3)-1)*r(3,:)-ref)
      table_pos(1,:)=dble(pos(1)-1)*r(1,:)+dble(pos(2)-1)*r(2,:)+dble(pos(3)-1)*r(3,:)-ref
      k=2

      do i=-1,1,1
       do j=-1,1,1
        table_dist(k)=norm(dble(pos(1)-1+i*dim_lat(1))*r(1,:)+dble(pos(2)-1+j*dim_lat(2))*r(2,:) &
         +dble(pos(3)-1)*r(3,:)-ref)
        table_pos(k,:)=dble(pos(1)+i*dim_lat(1)-1)*r(1,:)+dble(pos(2)+j*dim_lat(2)-1)*r(2,:)+ &
         dble(pos(3)-1)*r(3,:)-ref
        k=k+1
       enddo
      enddo

      min_dist=minloc(table_dist,1)
      pos_nei_2D=table_pos(min_dist,:)

      end function pos_nei_2D

!
! function that gives the position closest neighbor among the different supercell image
! with periodic boundary condition. The position might not be on a spin but between the spin
!
      function pos_nei_2D_real(tip,site,r,dim_lat)
      implicit none
      real(kind=8) :: pos_nei_2D_real(2)
      integer, intent(in) :: dim_lat(3)
      real(kind=8), intent(in) :: r(3,3),tip(2),site(2)
      !dummy
      integer :: i,j,k,min_dist
      real(kind=8) :: table_dist(10),table_pos(10,2),vec(2)

      vec=0.0d0

      k=1
      table_dist(1)=sqrt((tip(1)-site(1))**2+(tip(2)-site(2))**2)
      table_pos(1,:)=tip-site
      k=2

      do i=-1,1,1
       do j=-1,1,1
        vec=tip-dble(i*dim_lat(1))*r(1,1:2)-dble(j*dim_lat(2))*r(2,1:2)-site
        table_dist(k)=sqrt(vec(1)**2+vec(2)**2)
        table_pos(k,:)=vec(1:2)
        k=k+1
       enddo
      enddo

      min_dist=minloc(table_dist,1)
      pos_nei_2D_real=table_pos(min_dist,:)

      end function pos_nei_2D_real
!
! function that gives the position closest neighbor among the different supercell image
! with periodic boundary condition when the super-cell size is not known.
! The position might not be on a spin but between the spin
!
      function pos_nei_2D_user(tip,site,xmax,ymax)
      implicit none
      real(kind=8) :: pos_nei_2D_user(2)
      real(kind=8), intent(in) :: tip(2),site(2),xmax,ymax
      !dummy
      integer :: i,j,k,min_dist
      real(kind=8) :: table_dist(10),table_pos(10,2),vec(2)

      vec=0.0d0

      k=1
      table_dist(1)=sqrt((tip(1)-site(1))**2+(tip(2)-site(2))**2)
      table_pos(1,:)=tip-site
      k=2

      do i=-1,1,1
       do j=-1,1,1
        vec(1)=tip(1)-dble(i)*xmax-site(1)
        vec(2)=tip(2)-dble(j)*ymax-site(2)
        table_dist(k)=sqrt(vec(1)**2+vec(2)**2)
        table_pos(k,:)=vec(1:2)
        k=k+1
       enddo
      enddo

      min_dist=minloc(table_dist,1)
      pos_nei_2D_user=table_pos(min_dist,:)

      end function pos_nei_2D_user
!
! function that gives the position closest neighbor among the different supercell image
! with periodic boundary condition
!
      function pos_nei_3D_SL(ix,iy,iz,pos,r,dim_lat)
      use m_vector, only : norm
      implicit none
      real(kind=8) :: pos_nei_3D_SL(3)
      integer, intent(in) :: dim_lat(3),pos(3),ix,iy,iz
      real(kind=8), intent(in) :: r(3,3)
      !dummy
      integer :: i,j,k,min_dist
      real(kind=8) :: table_dist(10),table_pos(10,3),ref(3)

      ref=dble(ix-1)*r(1,:)+dble(iy-1)*r(2,:)+dble(iz-1)*r(3,:)
      k=1
      table_dist(1)=norm(dble(pos(1)-1)*r(1,:)+dble(pos(2)-1)*r(2,:)+dble(pos(3)-1)*r(3,:)-ref)
      table_pos(1,:)=dble(pos(1)-1)*r(1,:)+dble(pos(2)-1)*r(2,:)+dble(pos(3)-1)*r(3,:)-ref
      k=2

      do i=-1,1,1
       do j=-1,1,1
        table_dist(k)=norm(dble(pos(1)-1+i*dim_lat(1))*r(1,:)+dble(pos(2)-1+j*dim_lat(2))*r(2,:) &
         +dble(pos(3)-1)*r(3,:)-ref)
        table_pos(k,:)=dble(pos(1)+i*dim_lat(1)-1)*r(1,:)+dble(pos(2)+j*dim_lat(2)-1)*r(2,:)+ &
         dble(pos(3)-1)*r(3,:)-ref
        k=k+1
       enddo
      enddo

      min_dist=minloc(table_dist,1)
      pos_nei_3D_SL=table_pos(min_dist,:)

      end function pos_nei_3D_SL
!
! function that construct the rotation matrix of angle theta (in deg) with respect to arbitrary axis a
! source wikipedia.org/wiki/matrice_de_rotation
! +90 means that you are rotation anti clockwise
!
      function rot_mat_I(theta,a)
      use m_vector, only : norm
      use m_constants, only : pi
      implicit none
      real(kind=8) :: rot_mat_I(3,3)
      integer, intent(in) :: a(3)
      real(kind=8), intent(in) :: theta
      !dummy
      real(kind=8) :: c,s,an(3)

      c=cos(-theta*pi(1.0d0)/180.0d0)
      s=sin(-theta*pi(1.0d0)/180.0d0)
      an=dble(a)/norm(dble(a))

      rot_mat_I= &
      reshape((/an(1)**2+(1-an(1)**2)*c,an(1)*an(2)*(1-c)-an(3)*s,an(1)*an(3)*(1-c)+an(2)*s, &
                an(1)*an(2)*(1-c)+an(3)*s,an(2)**2+(1-an(2)**2)*c,an(2)*an(3)*(1-c)-an(1)*s, &
                an(1)*an(3)*(1-c)-an(2)*s,an(2)*an(3)*(1-c)-an(1)*s,an(3)**2+(1-an(3)**2)*c/) &
                ,(/3,3/))

      end function rot_mat_I

      function rot_mat_R(theta,a)
      use m_vector, only : norm
      use m_constants, only : pi
      implicit none
      real(kind=8) :: rot_mat_R(3,3),a(3)
      real(kind=8), intent(in) :: theta
      !dummy
      real(kind=8) :: c,s,an(3)

      c=cos(theta*pi(1.0d0)/180.0d0)
      s=sin(theta*pi(1.0d0)/180.0d0)
      an=a/norm(a)

      rot_mat_R= &
      reshape((/an(1)**2+(1-an(1)**2)*c,an(1)*an(2)*(1-c)-an(3)*s,an(1)*an(3)*(1-c)+an(2)*s, &
                an(1)*an(2)*(1-c)+an(3)*s,an(2)**2+(1-an(2)**2)*c,an(2)*an(3)*(1-c)-an(1)*s, &
                an(1)*an(3)*(1-c)-an(2)*s,an(2)*an(3)*(1-c)-an(1)*s,an(3)**2+(1-an(3)**2)*c/) &
                ,(/3,3/))

      end function rot_mat_R

      end module

! give all the possible corners of the direct space unit cell as a function of i_s
! the form of the output is
! ipu=is+corners(1)
! ipv=is+corners(2)
! ipuv=is+corners(3)
      subroutine givemecorners(r,zorder)
      use m_constants, only : pi
      use m_sym_utils
      use m_vector, only : norm
      implicit none
      integer  :: zorder
      real(kind=8), intent(in) :: r(3,3)
! dumy
      real(kind=8) :: rot_r(3,3),rotation(3,3),rot_ang,a(3),r_un(3,3),r_de(3,3)
      integer :: i,totcorn,au,av,k,avant

      do i=1,3
       a(i)=norm(r(i,:))
       r_un(i,:)=r(i,:)/a(i)
      enddo

!       corners(3*i+2,:)= corners(3*(i-1)+1,:)
!!! find rot_r(1) as a function of the r(1:2)
       if(dabs(r_un(1,2)*r_un(2,1)-r_un(1,1)*r_un(2,2)).lt.1.0d-7) then
        write(*,*) 'problem in findcorners'
        write(*,*) 'ay*bx-ax*by=0'
        stop
       endif

      if (mod(zorder,3).eq.0) then
       totcorn=6*3+3*3
      else
       totcorn=4*3
      endif
      allocate(corners(totcorn,3))

!               (/ipu,ipv,ipuv/) 

      i=0

      do while (i.lt.zorder)
       rot_ang=dble(i)*pi(2.0d0)/dble(zorder)
       rotation=reshape((/dcos(rot_ang),-dsin(rot_ang),0.0d0,&
        dsin(rot_ang),dcos(rot_ang),0.0d0,&
        0.0d0,0.0d0,1.0d0/),(/3,3/))

       rot_r=matmul(r_un,rotation)

       do k=1,2
        au=NINT(-(-rot_r(k,2)*r_un(2,1)+rot_r(k,1)*r_un(2,2))/(r_un(1,2)*r_un(2,1)-r_un(1,1)*r_un(2,2)))
        av=NINT(-(-rot_r(k,1)*r_un(1,2)+rot_r(k,2)*r_un(1,1))/(r_un(1,2)*r_un(2,1)-r_un(1,1)*r_un(2,2)))
        corners(3*i+k,:)= (/au,av,0/)
       enddo
       corners(3*i+3,:)=corners(3*i+1,:)+corners(3*i+2,:)
      i=i+1
      enddo

      if (mod(zorder,3).eq.0) then

       if (abs(zorder).eq.3) then
        r_de(1,:)=r_un(1,:)
        r_de(2,:)=r_un(1,:)+r_un(2,:)
        r_de(3,:)=r_un(3,:)    
        zorder=6
        avant=3*3
        else 
        r_de(1,:)=r_un(1,:)
        r_de(2,:)=r_un(1,:)-r_un(2,:)
        r_de(3,:)=r_un(3,:)
        zorder=3
        avant=6*3
       endif

       i=0

       do while (i.lt.abs(zorder))
        rot_ang=dble(i)*pi(2.0d0)/dble(zorder)
        rotation=reshape((/dcos(rot_ang),-dsin(rot_ang),0.0d0,&
         dsin(rot_ang),dcos(rot_ang),0.0d0,&
         0.0d0,0.0d0,1.0d0/),(/3,3/))

        rot_r=matmul(r_de,rotation)

        do k=1,2
         au=NINT(-(-rot_r(k,2)*r_un(2,1)+rot_r(k,1)*r_un(2,2))/(r_un(1,2)*r_un(2,1)-r_un(1,1)*r_un(2,2)))
         av=NINT(-(-rot_r(k,1)*r_un(1,2)+rot_r(k,2)*r_un(1,1))/(r_un(1,2)*r_un(2,1)-r_un(1,1)*r_un(2,2)))
         corners(avant+3*i+k,:)= (/au,av,0/)
        enddo
        corners(avant+3*i+3,:)=corners(avant+3*i+1,:)+corners(avant+3*i+2,:)
       i=i+1
       enddo
      endif

#ifdef CPP_DEBUG
      do i=1,totcorn
       write(*,*) corners(i,:)
      enddo
#endif
 
      end subroutine givemecorners
