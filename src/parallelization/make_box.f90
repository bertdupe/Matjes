      module m_make_box
      integer :: Xstart,Xstop,Ystart,Ystop,Zstart,Zstop
! size of the border for the MC with ghost
      logical,allocatable :: ghost_border(:,:,:,:)
! position of the up, down, right and left neighbor
      integer,allocatable :: rank_nei(:,:)
      end module m_make_box

      subroutine box_1D(Nx,n_ghost,Periodic_log,size_border)
      use m_make_box, only : Xstart,Xstop,Ystart,Ystop,Zstart,Zstop,ghost_border,rank_nei
      use m_mpi_prop
      implicit none
      integer, intent(in) :: Nx,n_ghost,size_border
      logical, intent(in) :: Periodic_log(3)
! internal variable
      integer :: reste,i,j,dims(2),periods,nei(4),l,max_rank,swap,ierr
      integer :: isize_box,MPI_COMM_WORKING
      integer, allocatable :: working_ranks(:),masters(:),trans(:)
      integer,allocatable :: box_ranks(:,:)
      logical :: discuss_trans,discuss

      include 'mpif.h'

! case 1D easier that the other, we do not have to care about the shape of the domain

      reste=1
      i=0
      Ystart=1
      Ystop=1
      Zstart=1
      Zstop=1
      discuss_trans=.False.
      discuss=.False.

      if (irank.eq.0) write(*,*) "1D domain decomposition"

      if ((n_ghost.ge.isize).and.(irank.eq.0)) then
       write(6,'(a)') "You must have more CPUs than domains"
       write(6,'(a)') "Increase the number of CPUS"
       call mpi_abort(MPI_COMM_WORLD,ierr)
      endif

      do while (reste.ne.0)
       reste=mod(Nx,n_ghost-i)
       if (i.ge.N_ghost) then
        write(6,'(a)') "problem in the domain decomposition"
        write(6,'(a)') "check make_box.f90"
        call mpi_abort(MPI_COMM_WORLD,ierr)
       endif
       i=i+1
      enddo

! update the number of proc which will actually work
! isize_box is the size of the domain decomposition.
! isize is the total number of procs: isize_box should divide isize. In case, it is not the case, this part is taking car of it.
      isize_box=N_ghost-i+1
      isize_working=isize/isize_box*isize_box
      allocate(box_ranks(isize/isize_box,isize_box))
      allocate(working_ranks(isize/isize_box*isize_box))

      do i=1,isize_working/isize_box
       do j=1,isize_box
        box_ranks(i,j)=isize_box*(i-1)+j-1
        working_ranks(isize_box*(i-1)+j)=isize_box*(i-1)+j-1
       enddo
      enddo
      max_rank=maxval(working_ranks)

! create a new group for the right amount of proc including the temperatures and the domains
      if (irank.le.max_rank) then
       call mpi_group_incl(all_world,isize_working,working_ranks,working_group,ierr)
       irank_working=irank
      else
       irank=-10
      endif
! ceate a new communicator for the entire working group (T+box)
      call mpi_comm_create(MPI_COMM_WORLD,working_group,MPI_COMM,ierr)
! create a new group for each boxes
      if (minval(abs(box_ranks(irank_working/isize_box+1,:)-irank_working)).eq.0) call mpi_group_incl(working_group,isize_box,box_ranks(irank_working/isize_box+1,:),working_box,ierr)
! create a new communicator for each working group
      call mpi_comm_create(MPI_COMM,working_box,MPI_COMM_BOX,ierr)
! define a new rank of the proc within each box
      call mpi_group_rank(working_box,irank_box,ierr)

      if (irank_working.eq.0) write(6,'(/,a,I6,a,/)') "we are actually using ", isize," processors"


! create the cartesian coordinate /temperature,Nx_box/
       dims=(/isize/isize_box,isize_box/)
! periodicity of the box
       if (Periodic_log(1)) then
        periods=1
       else
        periods=0
       endif

! create a cartesian coordinates of dimension dims with periodicity periods. The 0 correspond to no reordering the rank.
! I do not know exactly what it makes so I keep the ranks of irank_working
       call mpi_cart_create(MPI_COMM_BOX,1,isize_box,periods,0,MPI_COMM_CART,ierr)
! create the coordinates
       allocate(coords(1))

       call mpi_cart_coords(MPI_COMM_CART,irank_box,1,coords,ierr)

! attribute the starting coordinates
       Xstart=Nx/isize_box*coords(1)+1
       Xstop=Nx/isize_box*(coords(1)+1)

       isize=isize_working

! take care of the ghost domains
       N=(/Xstop-Xstart+1,Ystop-Ystart+1,Zstop-Zstart+1/)
       start=(/Xstart-1,Ystart-1,Zstart-1/)

! setup the border
       allocate(ghost_border(1,Xstart:Xstop,1,1))
       ghost_border=.False.

       ghost_border(1,Xstart:Xstart+size_border,1,1)=.True.
       ghost_border(1,Xstop-size_border:Xstop,1,1)=.True.

! store the neighbors and there position
       allocate(rank_nei(4,isize_working))
! order
! central rank=0
       nei(1)=irank_working
       nei(2)=irank_box
! direction 2 up=2, down=3
       call mpi_cart_shift(MPI_COMM_CART,0,1,nei(3),nei(4),ierr)
       call mpi_allgather(nei,4,MPI_INTEGER,rank_nei,4,MPI_INTEGER,MPI_COMM,ierr)

       write(6,'(2I4,a,I4,a,I4)') rank_nei(3:4,irank_working+1)," are neighbors of proc ", rank_nei(1,irank_working+1), " of box ", rank_nei(2,irank_working+1)
       call mpi_barrier(MPI_COMM,ierr)
       if (irank_working.eq.0) write(6,'(a)') "-----------------------"

! keeping tracks of all the masters
! the masters are the 0 rank for each box
       allocate(masters(isize/isize_box))
       allocate(trans(isize/isize_box))
       masters=0
       trans=0
       if (irank_box.eq.0) then
        trans(irank_working/isize_box+1)=irank_working
       endif
       call mpi_allreduce(trans,masters,isize/isize_box,MPI_INTEGER,MPI_SUM,MPI_COMM,ierr)

! create the group containing all masters
       call mpi_group_incl(working_group,isize/isize_box,masters,group_master,ierr)
! create the communicator associated with the master ranks of each replica
       call mpi_comm_create(MPI_COMM,group_master,MPI_COMM_MASTER,ierr)


       write(6,'(a,I6,a,I6,a,2I6)') "process ",irank_working," in the box ",irank_box," for spins in ",Xstart,Xstop

       call mpi_barrier(MPI_COMM,ierr)

#ifdef CPP_DEBUG
!!      call mpi_comm_rank(MPI_WORKING_WORLD,irank,ierr)
      write(6,*) irank,Xstart,Xstop
!      call mpi_finalize(ierr)
#endif

      end subroutine box_1D

!----------------------------------------------
      subroutine box_2D(Nx,Ny,n_ghost,Periodic_log,size_border)
      use m_make_box, only : Xstart,Xstop,Ystart,Ystop,Zstart,Zstop,ghost_border,rank_nei
      use m_mpi_prop
      use m_error
      implicit none
      integer, intent(in) :: Nx,Ny,n_ghost,size_border
      logical, intent(in) :: Periodic_log(3)
! internal variable
      integer :: reste,i,ierr,j,dims(2),periods(2),nei(6),errcode,k
      integer :: dim_proc(2),Nspin
      integer :: isize_box,MPI_COMM_WORKING,max_rank
      integer, allocatable :: working_ranks(:),masters(:),trans(:)
      integer,allocatable :: box_ranks(:,:),dim_proc_test(:,:)
      logical :: discuss_trans,discuss

      include 'mpif.h'

! case 1D easier that the other, we do not have to care about the shape of the domain

      reste=1
      i=0
      Zstart=1
      Zstop=1
      height=1

      if (irank.eq.0) write(*,*) 'entering 2D domain decomposition algo'

      if ((n_ghost.gt.isize).and.(irank.eq.0)) then
       write(6,'(a)') "You must have more CPUs than domains"
       write(6,'(a)') "Increase the number of CPUS"
       call mpi_abort(MPI_COMM_WORLD,errcode,ierr)
      endif

      do while (reste.ne.0)
       reste=mod(Nx*Ny,n_ghost-i)
       if (i.ge.isize) then
        write(6,'(a)') "problem in the domain decomposition"
        write(6,'(a)') "check make_box.f90"
        call mpi_abort(MPI_COMM_WORLD,errcode,ierr)
       endif
       i=i+1
      enddo

! update the number of proc which will actually work
      isize_box=n_ghost-i+1
      isize_working=isize/isize_box*isize_box
      allocate(box_ranks(isize/isize_box,isize_box))
      allocate(working_ranks(isize/isize_box*isize_box))
      allocate(dim_proc_test((NINT(sqrt(dble(isize_box)))+2)**2,2))

! number of proc per direction
      dim_proc=0
      dim_proc_test=0
      k=0
      do i=1,NINT(sqrt(dble(isize_box)))+2
       do j=1,NINT(sqrt(dble(isize_box)))+2
        if (i*j.eq.isize_box) then
         k=k+1
         dim_proc_test(k,:)=(/i,j/)
        endif
       enddo
      enddo

      if (k.eq.0) then
       write(6,'(a)') "chessboard domain decomposition failed"
       call MPI_ABORT(MPI_COMM_WORLD,ierr)
      endif

      i=dim_proc_test(k/2+1,1)
      j=dim_proc_test(k/2+1,2)

      if (i*j.ne.isize_box) then
       write(6,'(a)') "chessboard domain decomposition failed"
       call MPI_ABORT(MPI_COMM_WORLD,ierr)
      endif

      dim_proc=(/i,j/)
#ifdef CPP_DEBUG
      write(6,*) irank,dim_proc
#endif

! number of spin per domain
      Nspin=Nx*Ny/isize_box

! set the unused ranks to -10

      do i=1,isize_working/isize_box
       do j=1,isize_box
        box_ranks(i,j)=isize_box*(i-1)+j-1
        working_ranks(isize_box*(i-1)+j)=isize_box*(i-1)+j-1
       enddo
      enddo
      max_rank=maxval(working_ranks)

! create a new group for the right amount of proc including the temperatures and the domains
      if (irank.le.max_rank) then
       call mpi_group_incl(all_world,isize_working,working_ranks,working_group,ierr)
       irank_working=irank
      else
       irank=-10
      endif

! update the value of isize since isize=isize_working anyway
       isize=isize_working

! ceate a new communicator for the entire working group (T+box)
      call mpi_comm_create(MPI_COMM_WORLD,working_group,MPI_COMM,ierr)
! create a new group for each boxes
      if (minval(abs(box_ranks(irank_working/isize_box+1,:)-irank_working)).eq.0) call mpi_group_incl(working_group,isize_box,box_ranks(irank_working/isize_box+1,:),working_box,ierr)
! create a new communicator for each working group
      call mpi_comm_create(MPI_COMM,working_box,MPI_COMM_BOX,ierr)
      if (ierr.ne.0) call err_comm(ierr,MPI_COMM_BOX,'MPI_COMM_BOX','make_box')
! define a new rank of the proc within each box
      call mpi_group_rank(working_box,irank_box,ierr)

      if (irank_working.eq.0) write(6,'(/,a,I6,a,/)') "we are actually using ", isize," processors"

! calculate the domains. Look for something the more square possible

      do i=1,Nspin
       do j=-(i-1),(i-1)
        width=Nx/dim_proc(1)+j
        if (mod(Nspin,width).eq.0) exit
        if ((width.gt.Nx).or.(width.le.1)) then
         write(6,'(a)') "can not find domain decomposition"
         write(6,'(a)') "check box_2 or change number of procs"
         call MPI_ABORT(MPI_COMM_WORLD,ierr)
        endif
       enddo
       if (mod(Nspin,width).eq.0) exit
      enddo
! find the length

      do i=1,Nspin
       do j=-(i-1),(i-1)
        length=Ny/dim_proc(2)+j
        if (width*length.eq.Nspin) exit
        if (width*length.gt.Nspin) then
         write(6,'(a)') "can not find domain decomposition"
         write(6,'(a)') "check box_2 or change number of procs"
         call MPI_ABORT(MPI_COMM_WORLD,ierr)
        endif
       enddo
       if (width*length.eq.Nspin) exit
      enddo

      if ((((width.le.size_border).or.(length.le.size_border))).and.(irank.eq.0)) then
       write(6,'(a)') "the width or the length of the domain of the checkerboard decomposition"
       write(6,'(a)') "is smaller than the border"
       call MPI_ABORT(MPI_COMM_WORLD,ierr)
      endif

#ifdef CPP_DEBUG
      write(6,*) irank,"width=",width,"length=",length
#endif

! create the cartesian coordinate /temperature,Nx_box/
       dims=(/isize/isize_box,isize_box/)
! periodicity of the box
       if (Periodic_log(1)) then
        periods(1)=1
       else
        periods(1)=0
       endif
       if (Periodic_log(2)) then
        periods(2)=1
       else
        periods(2)=0
       endif

! create a cartesian coordinates of dimension dims with periodicity periods. The 1 correspond to the reordering the rank.
! the ranks are reodered in order to optimize the cummunications
       call mpi_cart_create(MPI_COMM_BOX,2,dim_proc,periods,1,MPI_COMM_CART,ierr)
! create the coordinates
       allocate(coords(2))

       call mpi_cart_coords(MPI_COMM_CART,irank_box,2,coords,ierr)

      Xstart=width*coords(1)+1
      Xstop=width*(coords(1)+1)
      if (Xstop.gt.Nx) Xstop=Nx

      Ystart=length*coords(2)+1
      Ystop=length*(coords(2)+1)
      if (Ystop.gt.Nx) Ystop=Ny

      if (irank_box.ge.isize_box) then
       Xstart=1
       Xstop=0
       Ystart=1
       Ystop=0
      endif

      write(6,'(2(a,I3),a,2I3)') "cordinates of proc ", irank_working," of domain ", irank_box," is ", coords

      if (irank.eq.0) write(6,'(a)') ""
! take care of the ghost domains
       N=(/Xstop-Xstart+1,Ystop-Ystart+1,Zstop-Zstart+1/)
       start=(/Xstart-1,Ystart-1,Zstart-1/)

! setup the border
       allocate(ghost_border(1,Xstart:Xstop,Ystart:Ystop,1))
       ghost_border=.False.

       ghost_border(1,Xstart:Xstart+size_border,:,1)=.True.
       ghost_border(1,Xstop-size_border:Xstop,:,1)=.True.
       ghost_border(1,:,Ystart:Ystart+size_border,1)=.True.
       ghost_border(1,:,Ystop-size_border:Ystop,1)=.True.

! store the neighbors and there position
       allocate(rank_nei(6,isize_working))
! order
! central rank=0
       nei(1)=irank_working
       nei(2)=irank_box
! The diagonal are not coded so far.
! direction 1 up=3, down=4
       call mpi_cart_shift(MPI_COMM_CART,0,1,nei(3),nei(4),ierr)
! direction 2 right=5, left=6
       call mpi_cart_shift(MPI_COMM_CART,1,1,nei(5),nei(6),ierr)
! diagonale
       call mpi_allgather(nei,6,MPI_INTEGER,rank_nei,6,MPI_INTEGER,MPI_COMM,ierr)

       call mpi_barrier(MPI_COMM,ierr)

       if (irank.eq.0) write(6,'(a)') "dimension x"
       write(6,'(2I4,a,I4,a,I4)') rank_nei(3:4,irank_working+1)," are neighbors of proc ", rank_nei(1,irank_working+1), " of domain ", rank_nei(2,irank_working+1)

       call mpi_barrier(MPI_COMM,ierr)

       if (irank.eq.0) write(6,'(a)') "dimension y"
       write(6,'(2I4,a,I4,a,I4)') rank_nei(5:6,irank_working+1)," are neighbors of proc ", rank_nei(1,irank_working+1), " of domain ", rank_nei(2,irank_working+1)
       call mpi_barrier(MPI_COMM,ierr)
       if (irank_working.eq.0) write(6,'(a)') "-----------------------"

! keeping tracks of all the masters
! the masters are the 0 rank for each box
       allocate(masters(isize/isize_box))
       allocate(trans(isize/isize_box))
       masters=0
       trans=0
       if (irank_box.eq.0) then
        trans(irank_working/isize_box+1)=irank_working
       endif
       call mpi_allreduce(trans,masters,isize/isize_box,MPI_INTEGER,MPI_SUM,MPI_COMM,ierr)

! create the group containing all masters
       call mpi_group_incl(working_group,isize/isize_box,masters,group_master,ierr)
! create the communicator associated with the master ranks of each replica
       call mpi_comm_create(MPI_COMM,group_master,MPI_COMM_MASTER,ierr)
       if (ierr.ne.0) call err_comm(ierr,MPI_COMM_MASTER,'MPI_COMM_MASTER','make_box')

       write(6,'(a,I6,a,I6,2(a,2I6))') "process ",irank_working," in the box ",irank_box, &
     &  " for spins along X ",Xstart,Xstop," and along Y ",Ystart,Ystop

       call mpi_barrier(MPI_COMM,ierr)

#ifdef CPP_DEBUG
      write(6,*) irank,"2D rank ","irank_nd","cartesian coords ", coords(1:2)
#endif

      end subroutine box_2D

!----------------------------------------------
      subroutine box_3D(Nx,Ny,Nz,n_ghost,Periodic_log)
      use m_make_box, only : Xstart,Xstop,Ystart,Ystop,Zstart,Zstop
      use m_mpi_prop
      implicit none
      integer, intent(in) :: Nx,Ny,Nz,n_ghost
      logical, intent(in) :: Periodic_log(3)
! internal variable
      integer :: reste,i,ierr,j,k,Nprocxy
      integer :: dim_proc(3),Nspin
      integer :: isize_box,MPI_COMM_WORKING
      integer,allocatable :: box_ranks(:,:)

      include 'mpif.h'

! case 1D easier that the other, we do not have to care about the shape of the domain

      reste=1
      i=0

      if (irank.eq.0) write(*,*) 'entering 3D domain decomposition algo'

      if (isize.lt.n_ghost) then
       write(6,'(a)') "You must have more CPUS than domains"
       write(6,'(a)') "Increase the number of CPUS"
       call mpi_abort(ierr)
      endif


      do while (reste.ne.0)
       reste=mod(Nx*Ny*Nz,n_ghost-i)
       if (i.ge.isize) then
        write(6,'(a)') "problem in the domain decomposition"
        write(6,'(a)') "check make_box.f90"
        call mpi_abort(MPI_COMM_WORLD,ierr)
       endif
       i=i+1
      enddo

! update the number of proc which will actually work
      isize_box=n_ghost-i+1
!      allocate(box_ranks(isize_box))
      do i=1,isize_box
!       box_ranks(i)=i-1
      enddo
#ifdef CPP_DEBUG
      if (irank.eq.0) write(6,*) box_ranks
#endif

! number of proc per direction
      dim_proc=0
      do i=1,NINT((dble(isize_box))**(1.0d0/3.0d0))+1
       do j=1,NINT((dble(isize_box))**(1.0d0/3.0d0))+1
        do k=1,NINT((dble(isize_box))**(1.0d0/3.0d0))
       if (i*j*k.eq.isize_box) exit
       if (i*j*k.gt.isize_box) then
        write(6,'(a)') "number of proc is a prime number"
        write(6,'(a)') "it should not be... How did you do that anyway!!"
        call MPI_ABORT(MPI_COMM_WORLD,ierr)
       endif
        enddo
        if (i*j*k.eq.isize_box) exit
       enddo
       if (i*j*k.eq.isize_box) exit
      enddo
      dim_proc=(/i,j,k/)
#ifdef CPP_DEBUG
      if (irank.eq.0) write(6,*) dim_proc
#endif

! number of spin per domain
      Nspin=Nx*Ny*Nz/isize_box
! create a new communicator for the right amount of proc
!      call mpi_group_incl(all_world,isize_box,box_ranks,working_group,ierr)
! create a new communicator for the group
      call mpi_comm_create(MPI_COMM_WORLD,working_group,MPI_COMM_WORKING,ierr)
! create a new cartesian communicator for the group
      call mpi_cart_create(MPI_COMM_WORKING,3,dim_proc,Periodic_log(1:3),.False.,MPI_COMM,ierr)

! calculate the domains. Look for something the more square possible

      do i=1,Nspin
       do j=-(i-1),(i-1)
        width=Nx/dim_proc(1)+j
        if (mod(Nspin,width).eq.0) exit
        if ((width.gt.Nx).or.(width.le.1)) then
         write(6,'(a)') "can not find domain decomposition along 1"
         write(6,'(a)') "check box_3 or change number of procs"
         call MPI_ABORT(MPI_COMM_WORLD,ierr)
        endif
       enddo
       if (mod(Nspin,width).eq.0) exit
      enddo
! find the length

      do i=1,Nspin
       do j=-(i-1),(i-1)
        length=Ny/dim_proc(2)+j
        if (mod(Nspin,width).eq.0) exit
        if ((length.gt.Nx).or.(length.le.1)) then
         write(6,'(a)') "can not find domain decomposition along 2"
         write(6,'(a)') "check box_3 or change number of procs"
         call MPI_ABORT(MPI_COMM_WORLD,ierr)
        endif
       enddo
       if (mod(Nspin,length).eq.0) exit
      enddo
! find the height

      do i=1,Nspin
       do j=-(i-1),(i-1)
        height=Nz/dim_proc(3)+j
        if (width*length*height.eq.Nspin) exit
        if (width*length*height.gt.Nspin) then
         write(6,'(a)') "can not find domain decomposition along 3"
         write(6,'(a)') "check box_3 or change number of procs"
         call MPI_ABORT(MPI_COMM_WORLD,ierr)
        endif
       enddo
       if (width*length*height.eq.Nspin) exit
      enddo
#ifdef CPP_DEBUG
      if (irank.eq.0) write(6,*) "width=",width,"length=",length,"height",height
#endif

      Xstart=width*mod(irank,dim_proc(1))+1
      Xstop=width*mod(irank+1,dim_proc(1))
      if ((mod(irank,dim_proc(1)).eq.(dim_proc(1)-1)).or.(dim_proc(1).eq.1)) Xstop=Nx

      Nprocxy=dim_proc(1)*dim_proc(2)
      Ystart=length*mod(mod(irank,Nprocxy)/dim_proc(1),dim_proc(2))+1
      Ystop=length*(mod(mod(irank,Nprocxy)/dim_proc(1),dim_proc(2))+1)
      if (dim_proc(2).eq.1) Ystop=Ny

      Zstart=height*(irank/Nprocxy)+1
      Zstop=height*(irank/Nprocxy+1)
      if (((irank/Nprocxy).eq.(dim_proc(3)-1)).or.(dim_proc(3).eq.1)) Zstop=Nz

      if (irank.ge.isize_box) then
       Xstart=1
       Xstop=0
       Ystart=1
       Ystop=0
       Zstart=1
       Zstop=0
      endif

      isize=isize_box

! setup the coordinates and associate a rank
      coords(1)=mod(irank,dim_proc(1))
      coords(2)=mod(mod(irank,Nprocxy)/dim_proc(1),dim_proc(2))
      coords(3)=irank/Nprocxy

!      call mpi_cart_rank(MPI_COMM,coords(1:3),irank_nd,ierr(1))

#ifdef CPP_DEBUG
      write(6,*) irank,"Xstart=",Xstart,"Xstop=",Xstop
      write(6,*) irank,"Ystart=",Ystart,"Ystop=",Ystop
      write(6,*) irank,"Zstart=",Zstart,"Zstop=",Zstop
#endif

      end subroutine box_3D
