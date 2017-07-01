       module m_setup_DM
       interface setup_DM
        module procedure setup_DM_1D
        module procedure setup_DM_2D
        module procedure setup_DM_3D
       end interface setup_DM
       contains

! routine that calculates the direction of the DM vector depending on the position of the non-magnetic atoms in the unit cell
! inputs
! ndm: number of DM vectors to calculate
! nei: number of shelves on which the DM are vectors are acting
! ind: index of the neighbors to consider
! r: lattice vectors
! motif: motif of atoms in the unit cell
! dim_lat: size of the supercell
! world: size of the world. The sze of the array matches the dimension of the problem: 1 for 1D, 2 for 2Ds...
! phase: 1 if thin films, 2 if multilayers
! c1: dummy to separate the cases
!
! output:
! a big array containing the DM vectors in cartesian coordinates


       function setup_DM_1D(ndm,nei,ind,r,motif,dim_lat,world,phase,c1)
       use m_vector, only : cross,norm
       use m_sym_utils, only : order_zaxis,angle_oriented,rot_mat
       use m_derived_types
       use m_constants, only : pi
       implicit none
! variable that come in
       integer, intent(in) :: c1
       integer, intent(in) :: dim_lat(3),nei,ndm,world(:),phase
       real (kind=8), intent(in) :: r(3,3)
       integer, intent(in) :: ind(nei)
       type (cell), intent(in) :: motif
! value of the function
       real (kind=8) :: setup_DM_1D(ndm,3,phase)
! dummy variable
       integer :: i,j,i_dm,n_z,k,step
       real(kind=8) :: vec(3),R1(3),R1_dum(3),R2(3),dumy,DM_vec(3),vec_2(3)
       real(kind=8) :: non_mag_at(count(.not.motif%i_m),3)
! part of the symmetry
       real(kind=8) :: rot_ang,rotation(3,3),angle
       real(kind=8) :: sign_z

       i_DM=0
       setup_DM_1D=0
       non_mag_at=0.0d0

! find none magnetic atoms in the motif
       i=0
       do j=1,size(motif%i_m)
        if (motif%i_m(j)) cycle
        i=i+1
        non_mag_at(i,:)=motif%pos(j,:)
       enddo
      ! let's order them from the closest to the (0,0,0) atom to the furthest.
       do j=1,i
        do k=1,i
         if (abs(non_mag_at(j,3)).gt.abs(non_mag_at(k,3))) then
         vec=non_mag_at(j,:)
         non_mag_at(j,:)=non_mag_at(k,:)
         non_mag_at(k,:)=vec
         endif
        enddo
       enddo

       if (count(motif%i_m).eq.1) then
        setup_DM_1D=0.0d0
       else
        do j=1,size(motif%i_m)
         if (.not.motif%i_m(j).and.(abs(sum(motif%pos(j,1:3))).lt.1.0d-8)) cycle
         do i=1,size(non_mag_at,1)
          R1=r(1,:)*non_mag_at(i,1)+r(2,:)*non_mag_at(i,2)+r(3,:)*non_mag_at(i,3)
          vec=r(1,:)*motif%pos(j,1)+r(2,:)*motif%pos(j,2)+r(3,:)*motif%pos(j,3)
          R2=R1+vec

          DM_vec=1/norm(R2)/norm(R1)/norm(R1-R2)*cross(R1,R2)
          setup_DM_1D(1,:,1)=DM_vec
          setup_DM_1D(2,2:3,1)=DM_vec(2:3)
          setup_DM_1D(2,1,1)=-DM_vec(1)
         enddo
        enddo
       endif

       end function setup_DM_1D

       function setup_DM_2D(ndm,nei,ind,r,motif,dim_lat,world,phase,c1,c2)
       use m_vector, only : cross,norm
       use m_sym_utils, only : order_zaxis,angle_oriented,rot_mat,pos_nei
       use m_derived_types
       use m_constants, only : pi
       use m_table_dist, only : Tdir
       implicit none
! variable that come in
       integer, intent(in) :: c1,c2
       integer, intent(in) :: dim_lat(3),nei,ndm,world(:),phase
       real (kind=8), intent(in) :: r(3,3)
       integer, intent(in) :: ind(nei)
       type (cell), intent(in) :: motif
! value of the function
       real (kind=8) :: setup_DM_2D(ndm,3,phase)
! dummy variable
       integer :: i,j,i_dm,n_z,k,step,i_nei,avant
       real(kind=8) :: vec(3),R1(3),R1_dum(3),R2(3),dumy,DM_vec(3),vec_2(3)
       real(kind=8) :: non_mag_at(count(.not.motif%i_m),3),mag_at(count(motif%i_m),3)
! directions of the neighbors
       real(kind=8) :: tabledir(3,nei)
! part of the symmetry
       real(kind=8) :: rot_ang,rotation(3,3),angle
       real(kind=8) :: sign_z
       logical :: exists

       step=1
       i_DM=0
       setup_DM_2D=0
       non_mag_at=0.0d0
       mag_at=0.0d0
       tabledir=0.0d0
       avant=0

! find none magnetic atoms in the motif
       i=0
       k=0
       do j=1,size(motif%i_m)
        if (motif%i_m(j))then
         k=k+1
         mag_at(k,:)=motif%pos(j,:)
        else
         i=i+1
         non_mag_at(i,:)=motif%pos(j,:)
        endif
       enddo

       call Tdir(r,nei,world,phase,motif,tabledir)

       n_z=order_zaxis(r)

        ! one magnetic atom in the unit cell, the rotation is done from u to v
       if ((count(motif%i_m).eq.1).or.(phase.ge.2)) then

         do i_nei=1,nei
          do i=1,size(non_mag_at,1)
           vec=r(1,:)*mag_at(1,1)+r(2,:)*mag_at(1,2)+r(3,:)*mag_at(1,3)
           vec_2=tabledir(:,i_nei)
           R1=r(1,:)*non_mag_at(i,1)+r(2,:)*non_mag_at(i,2)+r(3,:)*non_mag_at(i,3)
! find the position of the good non magnetic atom (goes to sym_utils.f90 for details)
           R1_dum=pos_nei(R1,vec,vec_2,r,dim_lat)

           R1=R1_dum-vec
           R2=R1_dum-vec_2
           DM_vec=cross(R1,R2)

            dumy=norm(DM_vec)
            DM_vec=DM_vec/dumy
             ! make sure that DM_vec and vec rotates in the sense of u toward v
!           write(*,*) 360.0d0/pi(2.0d0)*angle_oriented(vec,DM_vec)

           do j=1,ind(i_nei),2
            rot_ang=360.0d0/dble(n_z)*dble(j-1)
            rotation=rot_mat(rot_ang,(/0,0,1/))
            setup_DM_2D(avant+j,:,i)=matmul(rotation,DM_vec)
           enddo

           do j=2,ind(i_nei),2
            setup_DM_2D(avant+j,:,i)=-setup_DM_2D(avant+j-1,:,i)
           enddo

          enddo
          avant=avant+ind(i_nei)
         enddo

!!!!!!!!!!
!!!!!!!!!!
        elseif (count(motif%i_m).ne.1) then
         write(6,'(a)') "not coded"
        endif

       do j=1,phase
        do i=1,ndm
         dumy=norm(setup_DM_2D(i,:,j))
         if (dumy.gt.1.0d-8) setup_DM_2D(i,:,j)=setup_DM_2D(i,:,j)/dumy
        enddo
       enddo

       inquire(file='DM-2donly',exist=exists)
       if (exists) then
        setup_DM_2D(:,3,:)=0.0d0
        do j=1,phase
         do i=1,ndm
          dumy=norm(setup_DM_2D(i,:,j))
          setup_DM_2D(i,:,j)=setup_DM_2D(i,:,j)/dumy
         enddo
        enddo
       endif

#ifdef CPP_DEBUG
       write(*,*) ndm,nei,phase
       do j=1,phase
       i=1
       do while (i.le.ndm)
       write(*,*) setup_DM_2D(i,:,j)
       i=i+1
       enddo
       enddo
#endif

       end function setup_DM_2D

       function setup_DM_3D(ndm,nei,ind,r,motif,dim_lat,world,phase,c1,c2,c3)
       use m_vector, only : cross,norm
       use m_sym_utils, only : order_zaxis,angle_oriented,rot_mat
       use m_derived_types
       use m_constants, only : pi
       implicit none
! variable that come in
       integer, intent(in) :: c1,c2,c3
       integer, intent(in) :: dim_lat(3),nei,ndm,world(:),phase
       real (kind=8), intent(in) :: r(3,3)
       integer, intent(in) :: ind(nei)
       type (cell), intent(in) :: motif
! value of the function
       real (kind=8) :: setup_DM_3D(ndm,3,phase)
! dummy variable
       integer :: i,j,i_dm,n_z,k,step
       real(kind=8) :: vec(3),R1(3),R1_dum(3),R2(3),dumy,DM_vec(3),vec_2(3)
       real(kind=8) :: non_mag_at(count(.not.motif%i_m),3)
! part of the symmetry
       real(kind=8) :: rot_ang,rotation(3,3),angle
       real(kind=8) :: sign_z

       i_DM=0
       setup_DM_3D=0
       non_mag_at=0.0d0

! find none magnetic atoms in the motif
       i=0
       do j=1,size(motif%i_m)
        if (motif%i_m(j)) cycle
        i=i+1
        non_mag_at(i,:)=motif%pos(j,:)
       enddo
      ! let's order them from the closest to the (0,0,0) atom to the furthest.
       do j=1,i
        do k=1,i
         if (abs(non_mag_at(j,3)).gt.abs(non_mag_at(k,3))) then
         vec=non_mag_at(j,:)
         non_mag_at(j,:)=non_mag_at(k,:)
         non_mag_at(k,:)=vec
         endif
        enddo
       enddo

       if (count(motif%i_m).eq.1) then
        setup_DM_3D=0.0d0
       else
        do j=1,size(motif%i_m)
         if (.not.motif%i_m(j).and.(abs(sum(motif%pos(j,1:3))).lt.1.0d-8)) cycle
         do i=1,size(non_mag_at,1)
          R1=r(1,:)*non_mag_at(i,1)+r(2,:)*non_mag_at(i,2)+r(3,:)*non_mag_at(i,3)
          vec=r(1,:)*motif%pos(j,1)+r(2,:)*motif%pos(j,2)+r(3,:)*motif%pos(j,3)
          R2=R1+vec

          DM_vec=1/norm(R2)/norm(R1)/norm(R1-R2)*cross(R1,R2)
          setup_DM_3D(1,:,1)=DM_vec
          setup_DM_3D(2,2:3,1)=DM_vec(2:3)
          setup_DM_3D(2,1,1)=-DM_vec(1)
         enddo
        enddo
       endif

       end function setup_DM_3D

       end module
