      module m_reconstruct_mat
      interface rebuild_mat
       module procedure rebuild_2d
      end interface rebuild_mat
      contains

      subroutine rebuild_2d(spinafter,ncomm,spin)
      use m_mpi_prop, only : coords,width,length,height,irank,isize,MPI_COMM
      implicit none
      real(kind=8), intent(in) :: spinafter(:,:,:,:,:)
      real(kind=8), intent(inout) :: spin(:,:,:,:,:)
      integer, intent(in) :: ncomm
! internals
      real(kind=8) :: spintransfer(1:4,width,length,height,size(spinafter,5))
      integer :: local_coord(2),ierr(3)
      integer :: Xi,Xf,Yi,Yf,Zi,Zf
      integer :: i,i_x,i_y,i_z,i_m
      character(len=30) :: fname,toto

      include 'mpif.h'

      do i=0,isize-1

       spintransfer=spinafter
       local_coord=coords

! boradcast the spins
      call mpi_bcast(Spintransfer,ncomm,MPI_REAL8,i,MPI_COMM,ierr)
! broadcast the position of the spins
      call mpi_bcast(local_coord,3,MPI_INTEGER,i,MPI_COMM,ierr)

      Xi=width*local_coord(1)+1
      Xf=width*(local_coord(1)+1)
      Yi=length*local_coord(2)+1
      Yf=length*(local_coord(2)+1)
      Zi=1
      Zf=1

#ifdef CPP_DEBUG
      if (irank.eq.0) then
       write(6,*) i,"Xstart=",Xi,"Xstop=",Xf
       write(6,*) i,"Ystart=",Yi,"Ystop=",Yf
       write(6,*) i,"Zstart=",Zi,"Zstop=",Zf
      endif
#endif

#ifdef CPP_OPENMP
!$OMP parallel private(i_x,i_y,i_z,i_m) default(shared)
#endif
       do i_m=1,size(spinafter,5)
        do i_z=0,height-1
         do i_y=0,length-1
          do i_x=0,width-1
          spin(4:7,i_x+Xi,i_y+Yi,i_z+Zi,i_m)=Spintransfer(:,i_x+1,i_y+1,i_z+1,i_m)
          enddo
         enddo
        enddo
       enddo
#ifdef CPP_OPENMP
!$OMP end parallel
#endif
      enddo

#ifdef CPP_DEBUG
       write(fname,'(I3)') irank
       toto=trim(adjustl(fname))
       write(fname,'(a,18a,a)')'Spin_rebuild_',(toto(i:i),i=1, &
     & len_trim(toto)),'.dat'
       OPEN(15+irank,FILE=fname)
       do i_m=1,size(spinafter,5)
        do i_z=1,size(spin,4)
         do i_y=1,size(spin,3)
          do i_x=1,size(spin,2)
          write(15+irank,'(3f14.8)') spin(4:6,i_x,i_y,i_z,i_m)
          enddo
         enddo
        enddo
       enddo
       close(15+irank)
#endif
      end subroutine rebuild_2d

      end module m_reconstruct_mat
