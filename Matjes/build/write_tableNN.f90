       subroutine write_tableNN(tableNN,tot_N_Nneigh,dim_lat,n_m)
       implicit none
! input
       integer, intent(in) :: tot_N_Nneigh,dim_lat(3),n_m
       integer, intent(in) :: tableNN(4,tot_N_Nneigh,dim_lat(1),dim_lat(2),dim_lat(3),n_m)
! internal variables
       integer :: i,k,l,j,m


       open(7,file='write_tableNN.out',action='write',status='old',form='formatted')

       do i=1,n_m       ! number of atom in the unit cell (=1 in most cases)
         do j=1,dim_lat(3)      ! z coordinate
           do k=1,dim_lat(2)          !  y coordinate
             do l=1,dim_lat(1)              !  x coordinate
               do m=1,tot_N_Nneigh                 ! number of neighbor from 1 to tot_N_Nneigh

                 write(7,'(4I10)') tableNN(:,m,l,k,j,i)             ! 1: the neighbor has the same J    2: the neighbor is in the z direction

               enddo
             enddo
           enddo
         enddo
       enddo
       close(7)


       end subroutine write_tableNN
