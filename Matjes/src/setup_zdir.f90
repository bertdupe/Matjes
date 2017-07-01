      subroutine setup_zdir(phase,nvois,motif)
      use m_lattice, only : tableNN,indexNN,spin
      use m_rw_lattice, only : dim_lat,net
      use m_vector, only : norm
      use m_derived_types
#ifdef CPP_MPI
      use m_make_box, only : Xstart,Xstop,Ystart,Ystop,Zstart,Zstop
#endif
      implicit none
!!! inout var
      type (cell), intent(in) :: motif
      integer, intent(in) :: phase,nvois
!!! dummy
      integer :: j,i,k,l,i_m
!!! position of the 2 atoms
      real(kind=8) :: zs,zv
! position of the neighbors
      integer :: v_x,v_y,v_z,v_m
#ifndef CPP_MPI
      integer :: Xstart,Xstop,Ystart,Ystop,Zstart,Zstop
      Xstart=1
      Xstop=dim_lat(1)
      Ystart=1
      Ystop=dim_lat(2)
      Zstart=1
      Zstop=dim_lat(3)
#endif

!! super lattice case
      if (phase.eq.2) then
       do i_m=1,size(motif%i_m)
       if (.not.motif%i_m(i_m)) cycle
        do k=Zstart,Zstop
         do j=Ystart,Ystop
          do i=Xstart,Xstop

         l=1
          do while (l.le.nvois)
           v_x=tableNN(1,l,i,j,k,i_m)
           v_y=tableNN(2,l,i,j,k,i_m)
           v_z=tableNN(3,l,i,j,k,i_m)
           v_m=tableNN(4,l,i,j,k,i_m)

           zs=spin(3,i,j,k,i_m)
           zv=spin(3,v_x,v_y,v_z,v_m)
           if (dabs(zv-zs).gt.1.0d-7) then
            tableNN(5,l,i,j,k,i_m)=0
            tableNN(6,l,i,j,k,i_m)=1
           else
            tableNN(6,l,i,j,k,i_m)=0
            tableNN(5,l,i,j,k,i_m)=1
           endif
          l=l+1
          enddo
          enddo
         enddo
        enddo
       enddo
      else
       tableNN(5,:,:,:,:,:)=1
      endif

#ifdef CPP_DEBUG
      write(6,'(a)') 'in plane contribution'
      do i=1,dim_lat(1)
      do j=1,dim_lat(2)
      do k=1,dim_lat(3)
      do i_m=1,size(motif%i_m)
      if (motif%i_m(i_m)) cycle
       write(6,'('//repeat('x,I4',nvois)//')') (tableNN(5,l,i,j,k,i_m),l=1,nvois)
      enddo
      enddo
      enddo
      enddo

      write(6,'(a)') 'out of plane plane contribution'
      do i=1,dim_lat(1)
      do j=1,dim_lat(2)
      do k=1,dim_lat(3)
      do i_m=1,size(motif%i_m)
      if (motif%i_m(i_m)) cycle
       write(6,'('//repeat('x,I4',nvois)//')') (tableNN(6,l,i,j,k,i_m),l=1,nvois)
      enddo
      enddo
      enddo
      enddo
#endif

      end subroutine setup_zdir
