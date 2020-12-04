      subroutine order_lattice(spins,motif,my_lattice)
      use m_derived_types
      use m_vector, only : norm
      implicit none
      type(lattice), intent(in) :: my_lattice
      type(t_cell), intent(in) :: motif
      real(kind=8), intent(inout) :: spins(7,my_lattice%dim_lat(1),my_lattice%dim_lat(2),my_lattice%dim_lat(3),my_lattice%nmag)
! local variables
      real(kind=8) :: save_spin(7,my_lattice%dim_lat(1),my_lattice%dim_lat(2),my_lattice%dim_lat(3),my_lattice%nmag)
      integer :: i_x,i_y,i_z,n_x,n_y,n_z,dim_lat(3)
      real(kind=8) :: r_old(3,3),pos(3),min_dist(3),origin(3),denominator,r_u,r_v,r(3,3)

      r_old=0.0d0
      min_dist=0.0d0
      origin=0.0d0
      r=my_lattice%areal
      dim_lat=my_lattice%dim_lat
! first save the spin
      save_spin=spins

! find the old referential (must be cubic)
! we first find the motif in order to determine the index from the position
      pos=spins(1:3,1,1,1,1)

      min_dist(1)=minval(dabs(spins(1,:,:,:,:)-pos(1)),mask=dabs(spins(1,:,:,:,:)-pos(1)).gt.1.0d-8)
      min_dist(2)=minval(dabs(spins(2,:,:,:,:)-pos(2)),mask=dabs(spins(2,:,:,:,:)-pos(2)).gt.1.0d-8)
      if (dim_lat(3).ne.1) min_dist(3)=minval(dabs(spins(3,:,:,:,:)-pos(3)),mask=dabs(spins(3,:,:,:,:)-pos(3)).gt.1.0d-8)
! we find the origin of the old referential along x,y,z
      origin(1)=minval(dabs(spins(1,:,:,:,:)))
      origin(2)=minval(dabs(spins(2,:,:,:,:)))
      if (dim_lat(3).ne.1) origin(3)=minval(dabs(spins(3,:,:,:,:)))

! the old referential is then
      r_old(1,1)=min_dist(1)
      r_old(2,2)=min_dist(2)
      if (dim_lat(3).ne.1) then
       r_old(3,3)=min_dist(3)
       else
       r_old(3,3)=1.0d0
      endif

      if (dabs(r(1,2)*r(2,1)-r(1,1)*r(2,2)).lt.1.0d-8) then
       write(*,*) 'problem of lattice transfer in order_lattice.f90'
       stop
       else
       denominator=r(1,2)*r(2,1)-r(1,1)*r(2,2)
      endif

! calculate the new indexation
      do i_z=1,size(save_spin,4)
       do i_y=1,size(save_spin,3)
        do i_x=1,size(save_spin,2)
        ! old indices
!         o_x=nint((save_spin(i_x,i_y,i_z,1,1)-origin(1))/norm(r_old(1,:)))
!         o_y=nint((save_spin(i_x,i_y,i_z,1,2)-origin(2))/norm(r_old(2,:)))
!         o_x=mod(o_x-1+dim_lat(1),dim_lat(1))+1
!         o_y=mod(o_y-1+dim_lat(2),dim_lat(2))+1

         ! new indices
         r_u=(((save_spin(2,i_x,i_y,i_z,1)-origin(2))/norm(r_old(2,:))*r(2,1)- &
          (save_spin(1,i_x,i_y,i_z,1)-origin(1))/norm(r_old(1,:))*r(2,2))/denominator) &
          -motif%atomic(1)%position(1)+1

         r_v=((-(save_spin(2,i_x,i_y,i_z,1)-origin(2))/norm(r_old(2,:))*r(1,1)+ &
          (save_spin(1,i_x,i_y,i_z,1)-origin(1))/norm(r_old(1,:))*r(1,2))/denominator) &
          -motif%atomic(1)%position(2)+1


         n_z=i_z
         n_x=mod(nint(r_u)-1+dim_lat(1),dim_lat(1))+1
         n_y=mod(nint(r_v)-1+dim_lat(2),dim_lat(2))+1

         !update the positions
         save_spin(i_x,i_y,i_z,1,1:3)=r(1,:)*(dble(n_x-1)+motif%atomic(1)%position(1))+r(2,:)*(dble(n_y-1)+motif%atomic(1)%position(2))+ &
          r(3,:)*(dble(n_z-1)+motif%atomic(1)%position(3))

         ! transfer of the new spin lattice
         spins(n_x,n_y,n_z,:,:)=save_spin(i_x,i_y,i_z,:,:)
!         spins(n_x,n_y,n_z,:,4)=-save_spin(i_x,i_y,i_z,:,4)
!         spins(n_x,n_y,n_z,:,3)=-save_spin(i_x,i_y,i_z,:,3)

        enddo
       enddo
      enddo

      end subroutine order_lattice
