module m_arrange_neigh

private
public :: arrange_neigh

contains

subroutine arrange_neigh(DM_vector,tableNN,indexNN,dim_lat,net)
use m_sym_utils, only : pos_nei,rot_mat,order_zaxis
! The goal is to organize the neighbours in the same direction as the
! DM vectors
! the rotation starts with the u vector. The rotation goes then from u to v.
implicit none
! inout variables
real (kind=8), intent(in) :: net(3,3),DM_vector(:,:,:)
integer, intent(in) :: dim_lat(3),indexNN(:,:)
integer, intent(inout) :: tableNN(:,:,:,:,:,:)
! internal variable
! slopoes
integer :: ix,iy,iz,im,ig_x,ig_y,ig_z,ig_m
integer :: save_N(size(DM_vector,1),size(tableNN,1),size(tableNN,5))
! shapes
integer :: shape_DM_vectors(3),shape_index(2),shape_tableNN(6)
integer :: i,i_nei,k,avant,direct,test
real(kind=8) :: test_vec(3),position(3),test1,test2

avant=0
save_N=0
shape_DM_vectors=shape(DM_vector)
shape_tableNN=shape(tableNN)
shape_index=shape(indexNN)

direct=sign(1,order_zaxis(net))

do i_nei=1,shape_index(1)
  do im=1,shape_tableNN(6)
    do iz=1,shape_tableNN(5)
      do iy=1,shape_tableNN(4)
        do ix=1,shape_tableNN(3)


      do k=1,indexNN(i_nei,1)
        save_N(k,:,1)=tableNN(:,k+avant,ix,iy,iz,im)
      enddo

! find which neighbors is perpendicular to the first DM vector
      test=0
      do i=1,indexNN(i_nei,1)
        do k=1,indexNN(i_nei,1)

          position=pos_nei(ix,iy,iz,save_N(k,1:3,1),net,dim_lat)
          test_vec=matmul(rot_mat(-direct*90.0d0,(/0,0,1/)),DM_vector(avant+i,:,1))

          test1=dot_product(DM_vector(avant+i,:,1),position)
          test2=test_vec(1)*position(1)+test_vec(2)*position(2)

          if ((abs(test1).lt.1.0d-8).and.(test2.gt.0.0d0)) then
            tableNN(:,i+avant,ix,iy,iz,im)=save_N(k,:,1)
            test=test+1
            exit
          endif

        enddo
      enddo

      if (test.ne.indexNN(i_nei,1)) then
        write(6,'(a)') "problem with the DM"
        write(6,'(a,I2,a)') "test is ", test, "and should be 6"
        stop
      endif

        enddo
      enddo
    enddo
  enddo

avant=avant+indexNN(i_nei,1)

enddo


write(6,'(a)') ''
write(6,'(a)') 'DM vectors should match neighbors bounds'
write(6,*) "check for first site (1,1,1), (1,N/3,1), (N,N,1), (N/2,N/2,1), (N,N,N) only"

write(6,'(a)') '---------------------------'
write(6,'(a)') '(1,1,1)'
do k=1,size(DM_vector,1)
  write(6,*) DM_vector(k,:,1),tableNN(1:2,k,1,1,1,1)
enddo

if (dim_lat(2).ne.1) then
  write(6,'(a)') '---------------------------'
  write(6,'(a)') '(1,N/3,1)'
  do k=1,size(DM_vector,1)
    write(6,*) DM_vector(k,:,1),tableNN(1:2,k,1,dim_lat(2)/3,1,1)
  enddo
endif

write(6,'(a)') '---------------------------'
write(6,'(a)') '(N,N,1)'
do k=1,size(DM_vector,1)
  write(6,*) DM_vector(k,:,1),tableNN(1:2,k,dim_lat(1),dim_lat(2),1,1)
enddo

if (dim_lat(2).ne.1) then
  write(6,'(a)') '---------------------------'
  write(6,'(a)') '(N/2,N/2,1)'
  do k=1,size(DM_vector,1)
    write(6,*) DM_vector(k,:,1),tableNN(1:2,k,dim_lat(1)/2,dim_lat(2)/2,1,1)
  enddo
endif

if (dim_lat(3).ne.1) then
  write(6,'(a)') '---------------------------'
  write(6,'(a)') '(N/2,N/2,1)'
  do k=1,size(DM_vector,1)
    write(6,*) DM_vector(k,:,1),tableNN(1:2,k,dim_lat(1),dim_lat(2),dim_lat(3),1)
  enddo
endif

write(6,'(a)') ''


end subroutine arrange_neigh



!
!-------------------------------
!
end module m_arrange_neigh
