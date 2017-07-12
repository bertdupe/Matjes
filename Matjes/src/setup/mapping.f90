        module m_mapping
       interface mapping
        module procedure mapping_1D
        module procedure mapping_2D
        module procedure mapping_2D_SL
        module procedure mapping_2D_motif
        module procedure mapping_3D
        module procedure mapping_3D_motif_SL
       end interface mapping

       contains
! establish the table of neighbors
! it gives the indices of the neighbors of spin i_s
! mapping is a table which first contains
! first entry i_s: key of the cell
! second n is the neighbor index written on a line
! third column gives: 1-x;2-y;3-z position
       subroutine mapping_1D(r,n,d,nei,Nei_il,phase,motif,indexNN,tableNN)
       use m_rw_lattice, only : dim_lat
       use m_vector , only : norm
       use m_lattice, only : spin
       use m_derived_types
#ifdef CPP_MPI
       use m_make_box, only : Xstart
#endif
       implicit none
! variable that come in
       integer, intent(in) :: n,nei,Nei_il,phase
       type (cell), intent(in) :: motif
       real(kind=8), intent(in) :: d(:),r(3,3)
       integer, intent(in) :: indexNN(:)
! value of the function
       integer, intent(inout) :: tableNN(:,:,:)
! external blas
! 3D coordinate ix,iy,iz of 1d coordinate k
       integer :: i_x,i_y,i_z,i_m,Xstop
! dummy variable
       integer :: i_s,i,j,k,l,i_Nei,avant,i_p,i_phase
       real (kind=8) :: u,v,w,vec(3),dist
#ifndef CPP_MPI
       integer, parameter ::  Xstart=1
#endif

       avant=0
       Xstop=size(tableNN,3)

      do i_nei=1,nei
#ifdef CPP_OPENMP
!$OMP parallel private(i_x,l,i,i_p,vec,dist) default(shared)
#endif
        do i_x=0,Xstop-1

         l=1
         do i=-i_nei,i_nei,1
         do i_p=1,size(motif%i_m)
         if (.not.motif%i_m(i_p)) cycle
          vec=r(1,:)*(dble(i)+motif%pos(i_p,1))+r(2,:)*motif%pos(i_p,2)+ &
          r(3,:)*motif%pos(i_p,3)
          dist=norm(vec)

          if ((dabs(d(i_nei)-dist).lt.1.0d-8).and.(l.le.indexNN(i_Nei))) then
            tableNN(1,avant+l,i_x+1)=mod(i+i_x+Xstart+dim_lat(1)-1,dim_lat(1))+1
            tableNN(2,avant+l,i_x+1)=1
            tableNN(3,avant+l,i_x+1)=1
            tableNN(4,avant+l,i_x+1)=1
            l=l+1
           endif

          enddo
         enddo
        enddo
#ifdef CPP_OPENMP
!$OMP end parallel
#endif
        avant=avant+indexNN(i_Nei)
       enddo

#ifdef CPP_DEBUG
       do i_x=1,dim_lat(1)
         do l=1,size(tableNN,2)
          write(6,*) i_x,"i_x=",tableNN(1,l,i_x)
        enddo
       enddo
#endif

       end subroutine mapping_1D

!subroutine to treat the neighbour with one atom in the 2D unit cell
!
       subroutine mapping_2D(r,n,d,nei,Nei_il,phase,motif,indexNN,tableNN)
       use m_rw_lattice, only : dim_lat
       use m_vector , only : norm
       use m_lattice, only : spin
       use m_derived_types
#ifdef CPP_MPI
       use m_make_box, only : Xstart,Ystart
#endif
#ifdef CPP_OPENMP
       use omp_lib
#endif
       implicit none
       integer, intent(in) :: n,nei,Nei_il,phase
       type (cell), intent(in) :: motif
       real(kind=8), intent(in) :: d(:),r(3,3)
       integer, intent(in) :: indexNN(:)
! value of the function
       integer, intent(inout) :: tableNN(:,:,:,:)
! external blas
! 3D coordinate ix,iy,iz of 1d coordinate k
       integer :: i_x,i_y,i_z,i_m,Xstop,Ystop,Mstop
! size of the tableNN
       integer :: shape_tableNN(4)
! dummy variable
       integer :: i_s,i,j,k,l,i_Nei,avant,i_p,i_phase,transfer_nei
       real (kind=8) :: u,v,w,vec(3),dist
#ifndef CPP_MPI
       integer, parameter ::  Xstart=1
       integer, parameter ::  Ystart=1
#endif
#ifdef CPP_OPENMP
       integer :: ithread,nthreads

nthreads=omp_get_num_procs()
call omp_set_num_threads(nthreads)
#endif

       avant=0
       shape_tableNN=shape(tableNN)
       Xstop=shape_tableNN(3)
       Ystop=shape_tableNN(4)
       Mstop=size(motif%i_m)

#ifdef CPP_OPENMP
!!$OMP parallel default(shared) private(ithread)
!ithread=omp_get_thread_num()
#endif

        do i_nei=1,nei

         do i_y=0,Ystop-1
         do i_x=0,Xstop-1

           l=1
           do i=-i_nei,i_nei,1
           do j=-i_nei,i_nei,1
            do i_p=1,Mstop
           if (.not.motif%i_m(i_p)) cycle
           vec=r(1,:)*(dble(i)+motif%pos(i_p,1))+r(2,:)*(dble(j)+motif%pos(i_p,2))+ &
            r(3,:)*motif%pos(i_p,3)
           dist=norm(vec)

            if (dabs(d(i_nei)-dist).lt.1.0d-8) then
            tableNN(1,avant+l,i_x+1,i_y+1)=mod(i+i_x+Xstart+dim_lat(1)-1,dim_lat(1))+1
            tableNN(2,avant+l,i_x+1,i_y+1)=mod(j+i_y+Ystart+dim_lat(2)-1,dim_lat(2))+1
            tableNN(3,avant+l,i_x+1,i_y+1)=1
            tableNN(4,avant+l,i_x+1,i_y+1)=1
            l=l+1
            endif
            enddo
           enddo

          enddo
         enddo
        enddo
        avant=avant+indexNN(i_Nei)

       enddo

#ifdef CPP_OPENMP
!!$OMP end do
!!$OMP end parallel
#endif
#ifdef CPP_DEBUG
       do i_y=0,Ystop-1
        do i_x=0,Xstop-1
         do l=1,n
          write(6,*) i_x+Xstart,i_y+Ystart,"i_x=",tableNN(1,l,i_x+1,i_y+1),"i_y=",tableNN(2,l,i_x+1,i_y+1)
         enddo
        enddo
       enddo
#endif

       end subroutine mapping_2D

!subroutine to treat the SL case
!
       subroutine mapping_2D_SL(r,n,d,nei,Nei_il,phase,motif,indexNN,tableNN)
       use m_rw_lattice, only : dim_lat
       use m_vector , only : norm
       use m_lattice, only : spin
       use m_derived_types
#ifdef CPP_MPI
       use m_make_box, only : Xstart,Ystart
#endif
       implicit none
       integer, intent(in) :: n,nei,Nei_il,phase
       type (cell), intent(in) :: motif
       real(kind=8), intent(in) :: d(:,:),r(3,3)
       integer, intent(in) :: indexNN(:,:)
! value of the function
       integer, intent(inout) :: tableNN(:,:,:,:,:)
! external blas
! 3D coordinate ix,iy,iz of 1d coordinate k
       integer :: i_x,i_y,i_z,i_m,Xstop,Ystop
! dummy variable
       integer :: i_s,i,j,k,l,i_Nei,avant,i_p,i_phase
       real (kind=8) :: u,v,w,vec(3),dist
#ifndef CPP_MPI
       integer, parameter ::  Xstart=1
       integer, parameter ::  Ystart=1
#endif

       avant=0
       Ystop=size(tableNN,4)
       Xstop=size(tableNN,3)

        i_phase=1

        do i_nei=1,nei
#ifdef CPP_OPENMP
!$OMP parallel private(i_x,i_y,i_m,l,i,j,i_p,vec,dist) default(shared)
#endif
        do i_y=0,Ystop-1
         do i_x=0,Xstop-1
          do i_m=1,size(motif%i_m)
          if (.not.motif%i_m(i_m)) cycle

           l=1
           do i=-i_nei,i_nei,1
           do j=-i_nei,i_nei,1
           do i_p=1,size(motif%i_m)
           if (.not.motif%i_m(i_p)) cycle
           vec=r(1,:)*(dble(i)+motif%pos(i_p,1))+r(2,:)*(dble(j)+motif%pos(i_p,2))+ &
            r(3,:)*motif%pos(i_p,3)-(r(1,:)*motif%pos(i_m,1)+r(2,:)*motif%pos(i_m,2)+r(3,:)*motif%pos(i_m,3))
           dist=norm(vec)

           if (dabs(d(i_nei,i_phase)-dist).lt.1.0d-8) then
            tableNN(1,avant+l,i_x+1,i_y+1,i_m)=mod(i+i_x+Xstart+dim_lat(1)-1,dim_lat(1))+1
            tableNN(2,avant+l,i_x+1,i_y+1,i_m)=mod(j+i_y+Ystart+dim_lat(2)-1,dim_lat(2))+1
            tableNN(3,avant+l,i_x+1,i_y+1,i_m)=1
            tableNN(4,avant+l,i_x+1,i_y+1,i_m)=i_p
            l=l+1
           endif
           enddo
           enddo
           enddo

          enddo
         enddo
        enddo
#ifdef CPP_OPENMP
!$OMP end parallel
#endif
        avant=avant+indexNN(i_Nei,i_phase)
       enddo

        i_phase=2

        do i_nei=1,Nei_il
#ifdef CPP_OPENMP
!$OMP parallel private(i_x,i_y,i_m,l,i,j,i_p,vec,dist) default(shared)
#endif
        do i_y=0,Ystop-1
         do i_x=0,Xstop-1
          do i_m=1,size(motif%i_m)
          if (.not.motif%i_m(i_m)) cycle

           l=1
           do i=-i_nei,i_nei,1
           do j=-i_nei,i_nei,1
           do i_p=1,size(motif%i_m)
           if (.not.motif%i_m(i_p)) cycle
           vec=r(1,:)*(dble(i)+motif%pos(i_p,1))+r(2,:)*(dble(j)+motif%pos(i_p,2))+ &
            r(3,:)*motif%pos(i_p,3)-(r(1,:)*motif%pos(i_m,1)+r(2,:)*motif%pos(i_m,2)+r(3,:)*motif%pos(i_m,3))
           dist=norm(vec)

           if (dabs(d(i_nei,i_phase)-dist).lt.1.0d-8) then
            tableNN(1,avant+l,i_x+1,i_y+1,i_m)=mod(i+i_x+Xstart+dim_lat(1)-1,dim_lat(1))+1
            tableNN(2,avant+l,i_x+1,i_y+1,i_m)=mod(j+i_y+Ystart+dim_lat(2)-1,dim_lat(2))+1
            tableNN(3,avant+l,i_x+1,i_y+1,i_m)=1
            tableNN(4,avant+l,i_x+1,i_y+1,i_m)=i_p
            l=l+1
           endif
           enddo
           enddo
           enddo

          enddo
         enddo
        enddo
#ifdef CPP_OPENMP
!$OMP end parallel
#endif
        avant=avant+indexNN(i_Nei,i_phase)
       enddo

#ifdef CPP_DEBUG
       do i_x=1,Ystop
        do i_y=1,Xstop
         do i_m=1,size(motif%i_m)
         if (.not.motif%i_m(i_m)) cycle
         do l=1,size(tableNN,2)
          write(6,*) i_x,i_y,i_m,"i_x=",tableNN(1,l,i_x,i_y,i_m),"i_y=",tableNN(2,l,i_x,i_y,i_m), &
           "i_m=",tableNN(4,l,i_x,i_y,i_m)
         enddo
         enddo
        enddo
       enddo
#endif

       end subroutine mapping_2D_SL

!subroutine to treat more than one atom in the unit cell
!
       subroutine mapping_2D_motif(r,n,d,nei,Nei_il,phase,motif,indexNN,tableNN)
       use m_rw_lattice, only : dim_lat
       use m_vector , only : norm
       use m_lattice, only : spin
       use m_derived_types
#ifdef CPP_MPI
       use m_make_box, only : Xstart,Ystart
#endif
       implicit none
       integer, intent(in) :: n,nei,Nei_il,phase
       type (cell), intent(in) :: motif
       real(kind=8), intent(in) :: d(:),r(3,3)
       integer, intent(in) :: indexNN(:)
! value of the function
       integer, intent(inout) :: tableNN(:,:,:,:,:)
! external blas
! 3D coordinate ix,iy,iz of 1d coordinate k
       integer :: i_x,i_y,i_z,i_m,Xstop,Ystop
! dummy variable
       integer :: i_s,i,j,k,l,i_Nei,avant,i_p,i_phase
       real (kind=8) :: u,v,w,vec(3),dist
#ifndef CPP_MPI
       integer, parameter ::  Xstart=1
       integer, parameter ::  Ystart=1
#endif

       avant=0
       Xstop=size(tableNN,3)
       Ystop=size(tableNN,4)

        do i_nei=1,nei
#ifdef CPP_OPENMP
!$OMP parallel private(i_x,i_y,i_m,l,i,j,i_p,vec,dist) default(shared)
#endif
        do i_y=0,Ystop-1
         do i_x=0,Xstop-1
          do i_m=1,count(motif%i_m)

           l=1
           do i=-i_nei,i_nei,1
           do j=-i_nei,i_nei,1
           do i_p=1,size(motif%i_m)
           if (.not.motif%i_m(i_p)) cycle
           vec=r(1,:)*(dble(i)+motif%pos(i_p,1))+r(2,:)*(dble(j)+motif%pos(i_p,2))+ &
            r(3,:)*motif%pos(i_p,3)
           dist=norm(vec)

           if (dabs(d(i_nei)-dist).lt.1.0d-8) then
            tableNN(1,avant+l,i_x+1,i_y+1,i_m)=mod(i+i_x+Xstart+dim_lat(1)-1,dim_lat(1))+1
            tableNN(2,avant+l,i_x+1,i_y+1,i_m)=mod(j+i_y+Ystart+dim_lat(2)-1,dim_lat(2))+1
            tableNN(3,avant+l,i_x+1,i_y+1,i_m)=1
            tableNN(4,avant+l,i_x+1,i_y+1,i_m)=i_p
            l=l+1
           endif
           enddo
           enddo
           enddo

          enddo
         enddo
        enddo
#ifdef CPP_OPENMP
!$OMP end parallel
#endif
        avant=avant+indexNN(i_Nei)
       enddo

#ifdef CPP_DEBUG
       do i_x=1,dim_lat(1)
        do i_y=1,dim_lat(2)
         do i_m=1,size(motif%i_m)
         if (.not.motif%i_m(i_m)) cycle
         do l=1,size(tableNN,2)
          write(6,*) i_x,i_y,i_m,"i_x=",tableNN(1,l,i_x,i_y,i_m),"i_y=",tableNN(2,l,i_x,i_y,i_m), &
           "i_m=",tableNN(4,l,i_x,i_y,i_m)
         enddo
         enddo
        enddo
       enddo
#endif

       end subroutine mapping_2D_motif

       subroutine mapping_3D(r,n,d,nei,Nei_il,phase,motif,indexNN,tableNN)
       use m_rw_lattice, only : dim_lat
       use m_vector , only : norm
       use m_lattice, only : spin
       use m_derived_types
#ifdef CPP_MPI
       use m_make_box, only : Xstart,Ystart,Zstart
#endif
       implicit none
! variable that come in
       integer, intent(in) :: n,nei,Nei_il,phase
       type (cell), intent(in) :: motif
       real(kind=8), intent(in) :: d(:),r(3,3)
       integer, intent(in) :: indexNN(:)
! value of the function
       integer :: tableNN(:,:,:,:,:,:)
! external blas
! 3D coordinate ix,iy,iz of 1d coordinate k
       integer :: i_x,i_y,i_z,i_m,Xstop,Ystop,Zstop
! dummy variable
       integer :: i_s,i,j,k,l,i_Nei,avant,i_p,i_phase
       real (kind=8) :: u,v,w,vec(3),dist
#ifndef CPP_MPI
       integer, parameter ::  Xstart=1
       integer, parameter ::  Ystart=1
       integer, parameter ::  Zstart=1
#endif

       avant=0
       Xstop=size(tableNN,3)
       Ystop=size(tableNN,4)
       Zstop=size(tableNN,5)

       do i_nei=1,nei
#ifdef CPP_OPENMP
!$OMP parallel private(i_x,i_y,i_z,i_m,l,i,j,k,i_p,vec,dist) default(shared)
#endif
        do i_z=0,Zstop-1
         do i_y=0,Ystop-1
          do i_x=0,Xstop-1
           do i_m=1,count(motif%i_m)

! I selected a cell (i_x,i_y,i_z) and inside one of the atom i_m
! now we have to go along x,y and z direction and check all the distances
           l=1
           do i=-i_nei,i_nei,1
           do j=-i_nei,i_nei,1
           do k=-i_nei,i_nei,1
           do i_p=1,size(motif%i_m)
           if (.not.motif%i_m(i_p)) cycle

           vec=r(1,:)*(dble(i)+motif%pos(i_p,1))+r(2,:)*(dble(j)+motif%pos(i_p,2))+ &
            r(3,:)*(dble(k)+motif%pos(i_p,3))
           dist=norm(vec)

           if ((dabs(d(i_nei)-dist).lt.1.0d-8).and.(l.le.indexNN(i_Nei))) then
            tableNN(1,avant+l,i_x+1,i_y+1,i_z+1,i_m)=mod(i+i_x+Xstart+dim_lat(1)-1,dim_lat(1))+1
            tableNN(2,avant+l,i_x+1,i_y+1,i_z+1,i_m)=mod(j+i_y+Ystart+dim_lat(2)-1,dim_lat(2))+1
            tableNN(3,avant+l,i_x+1,i_y+1,i_z+1,i_m)=mod(k+i_z+Zstart+dim_lat(3)-1,dim_lat(3))+1
            tableNN(4,avant+l,i_x+1,i_y+1,i_z+1,i_m)=i_p
            l=l+1
           endif
           enddo
           enddo
           enddo
           enddo

           enddo
          enddo
         enddo
        enddo
#ifdef CPP_OPENMP
!$OMP end parallel
#endif
        avant=avant+indexNN(i_Nei)
       enddo


#ifdef CPP_DEBUG
       do i_x=1,dim_lat(1)
        do i_y=1,dim_lat(2)
         do i_z=1,dim_lat(3)
          do i_m=1,size(motif%i_m)
         if (.not.motif%i_m(i_m)) cycle
         do l=1,size(tableNN,2)
          write(6,*) i_x,i_y,i_z,i_m,"i_x=",tableNN(1,l,i_x,i_y,i_z,i_m),"i_y=",tableNN(2,l,i_x,i_y,i_z,i_m), &
           "i_z=",tableNN(3,l,i_x,i_y,i_z,i_m),"i_m=",tableNN(4,l,i_x,i_y,i_z,i_m)
         enddo
          enddo
         enddo
        enddo
       enddo
#endif

       end subroutine mapping_3D
!
!
!input variables are
!net,tot_N_Nneigh,tabledist(:,:),N_Nneigh,Nei_il,Nei_z,phase,motif,indexNN(:,:),tableNN(:,:,:,:,:,:)

       subroutine mapping_3D_motif_SL(r,n,d,nei,Nei_il,Nei_z,phase,motif,indexNN,tableNN)
       use m_rw_lattice, only : dim_lat
       use m_vector , only : norm
       use m_lattice, only : spin
       use m_derived_types
#ifdef CPP_MPI
       use m_make_box, only : Xstart,Ystart,Zstart
#endif
       implicit none
! variable that come in
       integer, intent(in) :: n,nei,Nei_il,phase,Nei_z
       type (cell), intent(in) :: motif
       real(kind=8), intent(in) :: d(:,:),r(3,3)
       integer, intent(in) :: indexNN(:,:)
! value of the function
       integer :: tableNN(:,:,:,:,:,:)
! external blas
! 3D coordinate ix,iy,iz of 1d coordinate k
       integer :: i_x,i_y,i_z,i_m,Xstop,Ystop,Zstop
! dummy variable
       integer :: i_s,i,j,k,l,i_Nei,avant,i_p,i_phase
       real (kind=8) :: u,v,w,vec(3),dist
#ifndef CPP_MPI
       integer, parameter ::  Xstart=1
       integer, parameter ::  Ystart=1
       integer, parameter ::  Zstart=1
#endif

       avant=0
       Xstop=size(tableNN,3)
       Ystop=size(tableNN,4)
       Zstop=size(tableNN,5)

        i_phase=1

        do i_nei=1,nei
#ifdef CPP_OPENMP
!$OMP parallel private(i_x,i_y,i_m,l,i,j,i_p,vec,dist) default(shared)
#endif
        do i_y=0,Ystop-1
         do i_x=0,Xstop-1
          do i_z=0,Zstop-1
           do i_m=1,size(motif%i_m)
           if (.not.motif%i_m(i_m)) cycle

            l=1
            do i=-i_nei,i_nei,1
            do j=-i_nei,i_nei,1
            do i_p=1,size(motif%i_m)
            if (.not.motif%i_m(i_p)) cycle
            vec=r(1,:)*(dble(i)+motif%pos(i_p,1))+r(2,:)*(dble(j)+motif%pos(i_p,2))+ &
             r(3,:)*motif%pos(i_p,3)-(r(1,:)*motif%pos(i_m,1)+r(2,:)*motif%pos(i_m,2)+r(3,:)*motif%pos(i_m,3))
            dist=norm(vec)

            if (dabs(d(i_nei,i_phase)-dist).lt.1.0d-8) then
             tableNN(1,avant+l,i_x+1,i_y+1,i_z+1,i_m)=mod(i+i_x+Xstart+dim_lat(1)-1,dim_lat(1))+1
             tableNN(2,avant+l,i_x+1,i_y+1,i_z+1,i_m)=mod(j+i_y+Ystart+dim_lat(2)-1,dim_lat(2))+1
             tableNN(3,avant+l,i_x+1,i_y+1,i_z+1,i_m)=i_z+Zstart
             tableNN(4,avant+l,i_x+1,i_y+1,i_z+1,i_m)=i_p
             l=l+1
            endif
            enddo
            enddo
            enddo

           enddo
          enddo
         enddo
        enddo
#ifdef CPP_OPENMP
!$OMP end parallel
#endif
        avant=avant+indexNN(i_Nei,i_phase)
       enddo

        i_phase=2

        do i_nei=1,Nei_il
#ifdef CPP_OPENMP
!$OMP parallel private(i_x,i_y,i_m,l,i,j,i_p,vec,dist) default(shared)
#endif
        do i_y=0,Ystop-1
         do i_x=0,Xstop-1
          do i_z=0,Zstop-1
           do i_m=1,size(motif%i_m)
           if (.not.motif%i_m(i_m)) cycle

            l=1
            do i=-i_nei,i_nei,1
            do j=-i_nei,i_nei,1
            do i_p=1,size(motif%i_m)
            if (.not.motif%i_m(i_p)) cycle
            vec=r(1,:)*(dble(i)+motif%pos(i_p,1))+r(2,:)*(dble(j)+motif%pos(i_p,2))+ &
             r(3,:)*motif%pos(i_p,3)-(r(1,:)*motif%pos(i_m,1)+r(2,:)*motif%pos(i_m,2)+r(3,:)*motif%pos(i_m,3))
            dist=norm(vec)

            if (dabs(d(i_nei,i_phase)-dist).lt.1.0d-8) then
             tableNN(1,avant+l,i_x+1,i_y+1,i_z+1,i_m)=mod(i+i_x+Xstart+dim_lat(1)-1,dim_lat(1))+1
             tableNN(2,avant+l,i_x+1,i_y+1,i_z+1,i_m)=mod(j+i_y+Ystart+dim_lat(2)-1,dim_lat(2))+1
             tableNN(3,avant+l,i_x+1,i_y+1,i_z+1,i_m)=i_z+Zstart
             tableNN(4,avant+l,i_x+1,i_y+1,i_z+1,i_m)=i_p
             l=l+1
            endif
            enddo
            enddo
            enddo

           enddo
          enddo
         enddo
        enddo
#ifdef CPP_OPENMP
!$OMP end parallel
#endif
        avant=avant+indexNN(i_Nei,i_phase)
       enddo

        i_phase=3

        do i_nei=1,Nei_z
#ifdef CPP_OPENMP
!$OMP parallel private(i_x,i_y,i_m,l,i,j,i_p,vec,dist) default(shared)
#endif
        do i_y=0,Ystop-1
         do i_x=0,Xstop-1
          do i_z=0,Zstop-1
           do i_m=1,size(motif%i_m)
           if (.not.motif%i_m(i_m)) cycle

            l=1
            do i=-i_nei,i_nei,1
            do j=-i_nei,i_nei,1
            do k=-i_nei,i_nei,1
            do i_p=1,size(motif%i_m)
            if (.not.motif%i_m(i_p)) cycle
            vec=r(1,:)*(dble(i)+motif%pos(i_p,1))+r(2,:)*(dble(j)+motif%pos(i_p,2))+ &
             r(3,:)*(dble(k)+motif%pos(i_p,3))-(r(1,:)*motif%pos(i_m,1)+r(2,:)*motif%pos(i_m,2)+r(3,:)*motif%pos(i_m,3))
            dist=norm(vec)

            if (dabs(d(i_nei,i_phase)-dist).lt.1.0d-8) then
             tableNN(1,avant+l,i_x+1,i_y+1,i_z+1,i_m)=mod(i+i_x+Xstart+dim_lat(1)-1,dim_lat(1))+1
             tableNN(2,avant+l,i_x+1,i_y+1,i_z+1,i_m)=mod(j+i_y+Ystart+dim_lat(2)-1,dim_lat(2))+1
             tableNN(3,avant+l,i_x+1,i_y+1,i_z+1,i_m)=mod(k+i_z+Zstart+dim_lat(3)-1,dim_lat(3))+1
             tableNN(4,avant+l,i_x+1,i_y+1,i_z+1,i_m)=i_p
             l=l+1
            endif
            enddo
            enddo
            enddo
            enddo

           enddo
          enddo
         enddo
        enddo
#ifdef CPP_OPENMP
!$OMP end parallel
#endif
        avant=avant+indexNN(i_Nei,i_phase)
       enddo


#ifdef CPP_DEBUG
       do i_x=1,dim_lat(1)
        do i_y=1,dim_lat(2)
         do i_z=1,dim_lat(3)
          do i_m=1,size(motif%i_m)
         if (.not.motif%i_m(i_m)) cycle
         do l=1,size(tableNN,2)
          write(6,*) i_x,i_y,i_z,i_m,"i_x=",tableNN(1,l,i_x,i_y,i_z,i_m),"i_y=",tableNN(2,l,i_x,i_y,i_z,i_m), &
           "i_z=",tableNN(3,l,i_x,i_y,i_z,i_m),"i_m=",tableNN(4,l,i_x,i_y,i_z,i_m),l

         enddo
          enddo
         enddo
        enddo
       enddo
#endif

       end subroutine mapping_3D_motif_SL

       end module m_mapping
