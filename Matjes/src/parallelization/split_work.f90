       module m_split_work
       interface split_work
        module procedure split_1D
        module procedure split_2D
        module procedure split_3D
       end interface split_work
       contains

! the only matrix that will not be split is the spin and the masque matrices. It simplifies a LOT the boundary issues.
! anyway what takes memory is the table of neighbors, the FFT setup and the predictors matrices.
! 1D
! create the different domains for the parallelisation of the code
       subroutine split_1D(dim_lat,Xstop,Xstart,isize,irank)
       use m_lattice, only : tableNN
       implicit none
       integer, intent(in) :: dim_lat(:),Xstop,Xstart,isize,irank
       !dumy
       integer, allocatable :: dumy(:,:,:,:,:,:)
       integer :: n_motif,n1,nx,ny,nz

       n_motif=size(tableNN,6)
       n1=size(tableNN,2)
       nx=dim_lat(1)
       ny=dim_lat(2)
       nz=dim_lat(3)

       allocate(dumy(6,n1,nx,ny,nz,n_motif))
       dumy=tableNN

! deallocate the table of neighbours and reallocate with the good size
       deallocate(tableNN)
       allocate(tableNN(6,n1,Xstart:Xstop,1,1,n_motif))

       tableNN=dumy(:,:,Xstart:Xstop,:,:,:)

       deallocate(dumy)

       end subroutine split_1D

! 2D
! create the different domains for the parallelisation of the code
       subroutine split_2D(dim_lat,Xstop,Xstart,Ystop,Ystart,isize,irank)
       use m_lattice, only : tableNN
       implicit none
       integer, intent(in) :: dim_lat(:),Xstop,Xstart,Ystop,Ystart,isize,irank
       !dumy
       integer, allocatable :: dumy(:,:,:,:,:,:)
       integer :: n_motif,n1,nx,ny,nz

       n_motif=size(tableNN,6)
       n1=size(tableNN,2)
       nx=dim_lat(1)
       ny=dim_lat(2)
       nz=dim_lat(3)

       allocate(dumy(6,n1,nx,ny,nz,n_motif))
       dumy=tableNN

! deallocate the table of neighbours and reallocate with the good size
       deallocate(tableNN)
       allocate(tableNN(6,n1,Xstart:Xstop,Ystart:Ystop,1,n_motif))

       tableNN=dumy(:,:,Xstart:Xstop,Ystart:Ystop,:,:)

       deallocate(dumy)

       end subroutine split_2D

! 3D
! create the different domains for the parallelisation of the code
       subroutine split_3D(dim_lat,Xstop,Xstart,Ystop,Ystart,Zstop,Zstart,isize,irank)
       use m_lattice, only : tableNN
       implicit none
       integer, intent(in) :: dim_lat(:),Xstop,Xstart,Ystop,Ystart,Zstop,Zstart,isize,irank
       !dumy
       integer, allocatable :: dumy(:,:,:,:,:,:)
       integer :: n_motif,n1,nx,ny,nz

       n_motif=size(tableNN,6)
       n1=size(tableNN,2)
       nx=dim_lat(1)
       ny=dim_lat(2)
       nz=dim_lat(3)

       allocate(dumy(6,n1,nx,ny,nz,n_motif))
       dumy=tableNN

! deallocate the table of neighbours and reallocate with the good size
       deallocate(tableNN)
       allocate(tableNN(6,n1,Xstart:Xstop,Ystart:Ystop,Zstart:Zstop,n_motif))

       tableNN=dumy(:,:,Xstart:Xstop,Ystart:Ystop,Zstart:Zstop,:)

       deallocate(dumy)

       end subroutine split_3D

       end module m_split_work
