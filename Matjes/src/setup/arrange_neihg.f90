       module m_arrange_neigh
       interface arrange_neigh
        module procedure arrange_2D
        module procedure arrange_2D_SL
        module procedure arrange_1D
        module procedure arrange_3D
        module procedure arrange_3D_SL
       end interface arrange_neigh

       contains

       subroutine arrange_2D(DM_vector,tableNN,indexNN,dim_lat,net,phase,tot_N_Nneigh,world,d,masque)
       use m_lattice, only : spin
       use m_parameters, only : DM
       use m_sym_utils, only : pos_nei,rot_mat,order_zaxis
#ifdef CPP_MPI
      use m_mpi_prop, only : start
#endif
       ! The goal is to organize the neighbours in the same direction as the
       ! DM vectors
       ! the rotation starts with the u vector. The rotation goes then from u to v.
       implicit none
       ! inout variables
       real (kind=8), intent(in) :: net(3,3),d(:),DM_vector(:,:)
       integer, intent(in) :: dim_lat(3),phase,tot_N_Nneigh,world(:),indexNN(:)
       integer, intent(inout) :: tableNN(:,:,:,:,:),masque(:,:,:)
       ! internal variable
       integer :: ix,iy,iz,save_N(size(DM_vector,1),size(tableNN,1),size(tableNN,5)),ig_x,ig_y
       integer :: save_masque(size(DM_vector,1))
       integer :: i,j,i_nei,k,avant,direct,Xstop,Ystop,test
       real(kind=8) :: test_vec(3),position(3),test1,test2
#ifndef CPP_MPI
       integer, dimension(3) :: start=0
#endif

       avant=0
       save_N=0
       direct=sign(1,order_zaxis(net))
       Xstop=size(tableNN,3)
       Ystop=size(tableNN,4)

       do i_nei=1,size(DM,1)
        do iy=1,Ystop
         do ix=1,Xstop
         ig_x=ix+start(1)
         ig_y=iy+start(2)
          do k=1,indexNN(i_nei)
          save_N(k,:,:)=tableNN(:,k+avant,ix,iy,:)
          save_masque(k)=masque(k+avant+1,ig_x,ig_y)
          enddo
! find which neighbors is perpendicular to the first DM vector
          test=0
          do i=1,indexNN(i_nei)
           do k=1,indexNN(i_nei)

            position=pos_nei(ig_x,ig_y,save_N(k,1:3,1),net,dim_lat)
            test_vec=matmul(rot_mat(-direct*90.0d0,(/0,0,1/)),DM_vector(avant+i,:))

            test1=dot_product(DM_vector(avant+i,:),position)
            test2=test_vec(1)*position(1)+test_vec(2)*position(2)

             if ((abs(test1).lt.1.0d-8).and.(test2.gt.0.0d0)) then
              tableNN(:,i+avant,ix,iy,:)=save_N(k,:,:)
              masque(i+avant+1,ig_x,ig_y)=save_masque(k)
              test=test+1
              exit
             endif

            enddo

          enddo
          if ((test.ne.6).and.(test.ne.4)) then
           write(6,'(a)') "problem with the DM"
           write(6,'(a,I2,a)') "test is ", test, "and should be 6"
           stop
          endif
         enddo
        enddo
        avant=avant+indexNN(i_nei)
       enddo

#ifdef CPP_DEBUG
       write(6,*) "DM vectors should match neighbors"
       write(6,*) "check for first site (1,1),(1,N/3),(N,N),(N/2,N/2)  only"
       do k=1,size(DM_vector,1)
        write(6,*) DM_vector(k,:),tableNN(1:2,k,1,1,1)
       enddo
       do k=1,size(DM_vector,1)
        write(6,*) DM_vector(k,:),tableNN(1:2,k,1,dim_lat(2)/3,1)
       enddo
       do k=1,size(DM_vector,1)
        write(6,*) DM_vector(k,:),tableNN(1:2,k,dim_lat(1),dim_lat(2),1)
       enddo
       do k=1,size(DM_vector,1)
        write(6,*) DM_vector(k,:),tableNN(1:2,k,dim_lat(1)/2,dim_lat(2)/2,1)
       enddo
#endif


       end subroutine arrange_2D
!
!-----------------------------------------------------------
!
       subroutine arrange_2D_SL(DM_vector,tableNN,indexNN,dim_lat,net,phase,tot_N_Nneigh,world,d,masque)
       use m_lattice, only : spin
       use m_parameters, only : DM
       use m_vector, only : norm
       use m_sym_utils, only : pos_nei,rot_mat
#ifdef CPP_MPI
       use m_mpi_prop, only : start
#endif
       ! The goal is to organize the neighbours in the same direction as the
       ! DM vectors
       ! the rotation starts with the u vector. The rotation goes then from u to v.
       implicit none
       ! inout variables
       real (kind=8), intent(in) :: net(3,3),d(:,:),DM_vector(:,:,:)
       integer, intent(in) :: dim_lat(3),phase,tot_N_Nneigh,world(:),indexNN(:,:)
       integer, intent(inout) :: tableNN(:,:,:,:,:),masque(:,:,:,:)
       ! internal variable
       integer :: ix,iy,iz,save_N(size(DM_vector,1),size(tableNN,1)),ig_x,ig_y,ig_z
       integer :: save_masque(size(DM_vector,1))
       integer :: i,j,i_nei,k,avant,Xstop,Ystop,i_m
       real(kind=8) :: test_vec(3),position(3),test1,test2
#ifndef CPP_MPI
       integer, dimension(3) :: start=0
#endif

       save_N=0
       Xstop=size(tableNN,3)
       Ystop=size(tableNN,4)

       do i_m=1,size(DM,2)
        avant=0
        do i_nei=1,size(DM(:,i_m),1)
         do iy=1,Ystop
          do ix=1,Xstop
           ig_x=ix+start(1)
           ig_y=iy+start(2)
           do k=1,indexNN(i_nei,1)
            save_N(k,:)=tableNN(:,k+avant,ix,iy,i_m)
            save_masque(k)=masque(k+avant+1,ig_x,ig_y,1)
           enddo

!          do k=1,indexNN(i_nei,1)
!           write(*,*) save_N(k,1:3),size(DM_vector,1)
!          enddo
!          pause
! find which neighbors is perpendicular to the first DM vector
           do i=1,indexNN(i_nei,1)
!           write(*,*) "------------------------------"
!           write(*,*) DM_vector(avant+i,:)

            do k=1,indexNN(i_nei,1)

             position=pos_nei(ig_x,ig_y,save_N(k,:),net,dim_lat)
             test_vec=matmul(rot_mat(-90.0d0,(/0,0,1/)),DM_vector(avant+i,:,i_m))

             test1=dot_product(DM_vector(avant+i,:,i_m),position)
             test2=(test_vec(1)*position(1)+test_vec(2)*position(2))/norm(position(1:2))/norm(test_vec(1:2))
!             write(*,*) ix,iy,position
!             write(*,*) save_N(k,:)
!             write(*,*) test1,test2
!             pause

             if ((abs(test1).lt.1.0d-8).and.(abs(test2-1.0d0).lt.1.0d-3)) then
              tableNN(:,avant+i,ix,iy,i_m)=save_N(k,:)
              masque(i+avant+1,ig_x,ig_y,1)=save_masque(k)
!              write(*,*) tableNN(1:4,avant+i,ix,iy,i_m),save_N(k,1:4)
!              write(*,*) DM_vector(avant+i,:,i_m),test_vec
!              write(*,*) i,k,ix,iy,i_m
!              pause
              exit
             endif

            enddo

           enddo
          enddo
         enddo
         avant=avant+indexNN(i_nei,1)
        enddo
       enddo

#ifdef CPP_DEBUG
       write(6,*) "DM vectors should match neighbors"
       write(6,*) "check for first site (1,1) only"
       write(6,*) size(DM_vector,1),size(DM_vector,3)
       do i=1,size(DM_vector,3)
       do k=1,size(DM_vector,1)
        write(6,*) DM_vector(k,:,i),tableNN(1:4,k,1,1,i)
       enddo
       enddo
       stop
#endif
       end subroutine arrange_2D_SL
!
!-------------------------------
!
! variable comming in
! DM_vector,tableNN(:,:,:,:,:,:),indexNN(:,:),dim_lat,net,phase,tot_N_Nneigh,world,tabledist(:,:),masque

       subroutine arrange_3D_SL(DM_vector,tableNN,indexNN,dim_lat,net,phase,tot_N_Nneigh,world,d,masque)
       use m_lattice, only : spin
       use m_parameters, only : DM
       use m_sym_utils, only : pos_nei,rot_mat
#ifdef CPP_MPI
       use m_mpi_prop, only : start
#endif
!       ! The goal is to organize the neighbours in the same direction as the
       ! DM vectors
       ! the rotation starts with the u vector. The rotation goes then from u to v.
       implicit none
       ! inout variables
       real (kind=8), intent(in) :: net(3,3),d(:,:),DM_vector(:,:)
       integer, intent(in) :: dim_lat(3),phase,tot_N_Nneigh,world(:),indexNN(:,:)
       integer, intent(inout) :: tableNN(:,:,:,:,:,:),masque(:,:,:,:)
       ! internal variable
       integer :: ix,iy,iz,save_N(size(DM_vector,1),size(tableNN,1),size(tableNN,6)),ig_x,ig_y,ig_z
       integer :: save_masque(size(DM_vector,1))
       integer :: i,j,i_nei,k,avant,Xstop,Ystop,Zstop
       real(kind=8) :: test_vec(3),position(3),test1,test2
#ifndef CPP_MPI
       integer, dimension(3) :: start=0
#endif

       avant=0
       save_N=0
       Xstop=size(tableNN,3)
       Ystop=size(tableNN,4)
       Zstop=size(tableNN,5)

       do i_nei=1,size(DM)
        do iy=1,Ystop
         do ix=1,Xstop
          do iz=1,Zstop
           ig_x=ix+start(1)
           ig_y=iy+start(2)
           ig_z=iz+start(3)

         do k=1,indexNN(i_nei,1)
         save_N(k,:,:)=tableNN(:,k+avant,ix,iy,iz,:)
         save_masque(k)=masque(k+avant+1,ig_x,ig_y,ig_z)
         enddo

!         do k=1,indexNN(i_nei,1)
!          write(*,*) save_N(k,1,1:3),size(DM_vector,1)
!         enddo
!         pause
! find which neighbors is perpendicular to the first DM vector
          do i=1,indexNN(i_nei,1)
!          write(*,*) "------------------------------"
!          write(*,*) DM_vector(avant+i,:)

            do k=1,indexNN(i_nei,1)

             position=pos_nei(ig_x,ig_y,ig_z,save_N(k,:,1),net,dim_lat)
             test_vec=matmul(rot_mat(-90.0d0,(/0,0,1/)),DM_vector(avant+i,:))

             test1=dot_product(DM_vector(avant+i,:),position)
             test2=test_vec(1)*position(1)+test_vec(2)*position(2)

             if ((dabs(test1).lt.1.0d-8).and.(dabs(test2-1.0d0).lt.1.0d-3)) then
              tableNN(:,avant+i,ix,iy,iz,:)=save_N(k,:,:)
              masque(i+avant+1,ig_x,ig_y,ig_z)=save_masque(k)
              exit
             endif

            enddo
           enddo
          enddo
         enddo
        enddo
        avant=avant+indexNN(i_nei,1)
       enddo

#ifdef CPP_DEBUG
       write(6,*) "DM vectors should match neighbors"
       write(6,*) "check for first site (1,1) only"
       do k=1,size(DM_vector,1)
        write(6,*) DM_vector(k,:),tableNN(1:2,k,1,1,1,1)
       enddo
#endif

       end subroutine arrange_3D_SL
!
!-------------------------------
!
       subroutine arrange_3D(DM_vector,tableNN,indexNN,dim_lat,net,phase,tot_N_Nneigh,world,d,masque)
       use m_lattice, only : spin
       ! The goal is to organize the neighbours in the same direction as the
       ! DM vectors
       ! the rotation starts with the u vector. The rotation goes then from u to v.
       implicit none
       ! inout variables
       real (kind=8), intent(in) :: net(3,3),d(:),DM_vector(:,:)
       integer, intent(in) :: dim_lat(3),phase,tot_N_Nneigh,world(:),indexNN(:)
       integer, intent(inout) :: tableNN(:,:,:,:,:,:),masque(:,:,:,:)
       ! internal variable
       integer :: ix,iy,iz
       integer :: i,j
       real(kind=8) :: test(3)

       write(*,*) "arrange neigh not coded"
       stop
       end subroutine arrange_3D

!
!-------------------------------
!
       subroutine arrange_1D(DM_vector,tableNN,indexNN,dim_lat,net,phase,tot_N_Nneigh,world,d,masque)
       use m_lattice, only : spin
       ! The goal is to organize the neighbours in the same direction as the
       ! DM vectors
       ! the rotation starts with the u vector. The rotation goes then from u to v.
       implicit none
       ! inout variables
       real (kind=8), intent(in) :: net(3,3),d(:),DM_vector(:,:)
       integer, intent(in) :: dim_lat(3),phase,tot_N_Nneigh,world(:),indexNN(:)
       integer, intent(inout) :: tableNN(:,:,:,:),masque(:,:,:,:)
       ! internal variable
       integer :: ix,iy,iz
       integer :: i,j
       real(kind=8) :: test(3)

       write(*,*) "arrange neigh not coded"
       stop
       end subroutine arrange_1D

       end module m_arrange_neigh
