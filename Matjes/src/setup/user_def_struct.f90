      subroutine user_def_struct(masque,spins,dim_lat,N_M,Nnei,net)
      use m_vector, only : norm,AXeqB
      use m_lattice, only : tableNN
#ifdef CPP_MPI
      use m_mpi_prop, only : irank
      use m_make_box, only : Xstart,Xstop,Ystart,Ystop,Zstart,Zstop
#endif
      implicit none
      integer, intent(in) :: dim_lat(3),N_M,Nnei
      real(kind=8), intent(in) :: net(3,3)
      real(kind=8), intent(inout) :: spins(7,dim_lat(1),dim_lat(2),dim_lat(3),N_M)
#ifdef CPP_MPI
      integer, intent(inout) :: masque(Nnei,Xstart:Xstop,Ystart:Ystop,Zstart:Zstop)
#else
      integer, intent(inout) :: masque(Nnei,dim_lat(1),dim_lat(2),dim_lat(3))
#endif
! internal variable
      integer, parameter :: io=11
      character(len=2) :: dumy
! position of the atom in the pattern
      real(kind=8), allocatable :: position(:,:)
! coordinate of the atom in the pattern
      integer, allocatable :: coord(:,:)
! number of atoms in the file
      integer :: natom
! number of cell and spin in the spin file
      integer :: ncell,nspin
! stuff
      integer :: i,j,i_x,i_y,i_z,i_m,i_pos,vecori(3)
      integer :: v_x,v_y,v_z
      real(kind=8) :: origine(3),lattice(3,3),r(3,3)
      integer :: upper(3),lower(3)

! put all the onsite masque value to 0
#ifdef CPP_MPI
      do i_z=Zstart,Zstop
       do i_y=Ystart,Ystop
        do i_x=Xstart,Xstop
         masque(1,i_x,i_y,i_z)=0
        enddo
       enddo
      enddo
#else
      do i_z=1,dim_lat(3)
       do i_y=1,dim_lat(2)
        do i_x=1,dim_lat(1)
         masque(1,i_x,i_y,i_z)=0
        enddo
       enddo
      enddo
#endif
#ifdef CPP_MPI
      if (irank.eq.0) then
      write(6,'(a)') 'user defined structure is read'
      write(6,'(a)') 'the 2 first lines of structure.xyz MUST BE the good lattice!!'
      endif
#else
      write(6,'(a)') 'user defined structure is read'
      write(6,'(a)') 'the 2 first lines of structure.xyz MUST BE the good lattice!!'
#endif

!calculate the number of spin and cell
      ncell=product(dim_lat)
      nspin=ncell*N_M

! lattice of lattice.in
      do i=1,3
       r(i,:)=net(i,:)/norm(net(i,:))
      enddo

      open(io,file='structure.xyz',action='read',status='old',form='formatted')
      rewind(io)
      do i=1,3
       read(io,*) (lattice(i,j),j=1,3)
      enddo
      read(io,*) natom
      allocate(position(3,natom))
      read(io,*) dumy
      do i=1,natom
       read(io,*) dumy,(position(j,i),j=1,3)
      enddo
      close(io)

! check that the number of cell is greater that the number of atoms in the file
      if (ncell.lt.natom) then
       write(*,*) 'the super cell set up in lattice.in is too small'
       stop
      endif

! the coordinate of the atoms are put into a lattice
! be carefull, the all start at 0,0,0 of course
      origine=position(:,1)
      vecori=AXeqB(transpose(lattice),position(:,1))
      allocate(coord(3,natom))
      do i_pos=1,natom
       coord(:,i_pos)=AXeqB(transpose(lattice),position(:,i_pos))
       if (norm(position(:,i_pos)).lt.norm(origine)) then
        origine=position(:,i_pos)
        vecori=coord(:,i_pos)
       endif
      enddo
! put the shape into a box
      do i=1,3
       upper(i)=maxval(coord(i,:))
       lower(i)=minval(coord(i,:))
      enddo
! update the index
      do i_pos=1,natom
       coord(:,i_pos)=coord(:,i_pos)-lower
      enddo
      upper=upper-lower
      lower=lower-lower


#ifdef CPP_DEBUG
      write(*,*) upper
      write(*,*) lower
#endif

! now we can update the onsite masque safely safely
      do i_pos=1,natom
       i_x=coord(1,i_pos)+1
       i_y=coord(2,i_pos)+1
       i_z=coord(3,i_pos)+1
       masque(1,i_x,i_y,i_z)=1
!#ifdef CPP_DEBUG
       write(*,*) i_x,i_y,i_z,masque(1,i_x,i_y,i_z)
!#endif
      enddo
! we have to update the spins and the masque of the neighbours to take into account the boundaries
! and defects

#ifdef CPP_MPI
      do i_z=Zstart,Zstop
       do i_y=Ystart,Ystop
        do i_x=Xstart,Xstop
#else
      do i_z=1,dim_lat(3)
       do i_y=1,dim_lat(2)
        do i_x=1,dim_lat(1)
#endif
         do i=2,Nnei
         ! if no atom is on site i_x,i_y,i_z then its interaction with the neighbours
         ! must be zero
          if (masque(1,i_x,i_y,i_z).eq.0) then
           spins(4:7,i_x,i_y,i_z,:)=0.0d0
           masque(:,i_x,i_y,i_z)=0
          else
         ! if there is an atom on site i_x,i_y,i_z we have to check if there is an atom
         ! on the neighbours sites
           do i_m=1,N_M
           v_x=tableNN(1,i-1,i_x,i_y,i_z,i_m)
           v_y=tableNN(2,i-1,i_x,i_y,i_z,i_m)
           v_z=tableNN(3,i-1,i_x,i_y,i_z,i_m)
            if (masque(1,v_x,v_y,v_z).eq.0) masque(i,i_x,i_y,i_z)=0
           enddo
          endif

         enddo
        enddo
       enddo
      enddo

#ifdef CPP_DEBUG
      do i_z=1,dim_lat(3)
       do i_y=1,dim_lat(2)
        do i_x=1,dim_lat(1)
       write(*,*) i_x,i_y,i_z,masque(:,i_x,i_y,i_z)
        enddo
       enddo
      enddo
#endif

      deallocate(position)
      end subroutine user_def_struct
