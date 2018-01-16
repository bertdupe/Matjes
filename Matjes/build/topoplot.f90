      module m_topoplot
      interface topoplot
       module procedure topoplot_2D
       module procedure topoplot_3D
       module procedure topoplot_1D
      end interface
      contains

      subroutine topoplot_1D(signature,map_eul,map_vort)

      implicit none
      real(kind=8), intent(in) :: map_eul(:)
      real(kind=8),intent(in) :: map_vort(:,:)
      real(kind=8), intent(in) :: signature
! inteernal variables
      integer :: i_x,j,i
      integer :: shape_map(1)
      character(len=30) :: fname,toto

      shape_map=shape(map_eul)



      write(6,'(a)') 'plot the topological charge'

      write(fname,'(f8.4)') signature
      toto=trim(adjustl(fname))
      write(fname,'(a,18a,a)')'topomap_',(toto(i:i),i=1, &
       len_trim(toto)),'.dat'
      OPEN(70,FILE=fname,action='write',status='unknown',form='formatted')

      do i_x=1,shape_map(1)

        Write(70,'((I6,2x),4(f14.8,2x))') i_x,(map_vort(j,i_x),j=1,3), &
     & map_eul(i_x)

      enddo
      close(70)

      end subroutine topoplot_1D

      subroutine topoplot_2D(signature,map_eul,map_vort)
      implicit none
      real(kind=8), intent(in) :: map_eul(:,:)
      real(kind=8),intent(in) :: map_vort(:,:,:)
      real(kind=8), intent(in) :: signature
! inteernal variables
      integer :: i_x,i_y,j,i
      integer :: shape_map(2)
      character(len=30) :: fname,toto

      shape_map=shape(map_eul)

      write(6,'(a)') 'plot the topological charge'

      write(fname,'(f8.4)') signature
      toto=trim(adjustl(fname))
      write(fname,'(a,18a,a)')'topomap_',(toto(i:i),i=1, &
       len_trim(toto)),'.dat'
      OPEN(70,FILE=fname,action='write',status='unknown',form='formatted')

      do i_y=1,shape_map(2)
       do i_x=1,shape_map(1)

        Write(70,'(2(I6,2x),4(f14.8,2x))') i_x,i_y,(map_vort(j,i_x,i_y),j=1,3), &
     & map_eul(i_x,i_y)

       enddo
      enddo
      close(70)

      end subroutine topoplot_2D

      subroutine topoplot_3D(signature,map_vort,map_eul)

      implicit none
      real(kind=8), intent(in) :: map_eul(:,:,:)
      real(kind=8),intent(in) :: map_vort(:,:,:,:)
      real(kind=8), intent(in) :: signature
! inteernal variables
      integer :: i_x,i_y,i_z,j,i
      integer :: shape_map(3)
      character(len=30) :: fname,toto

      shape_map=shape(map_eul)

      write(6,'(a)') 'plot the topological charge'

      write(fname,'(f8.4)') signature
      toto=trim(adjustl(fname))
      write(fname,'(a,18a,a)')'topomap_',(toto(i:i),i=1, &
       len_trim(toto)),'.dat'
      OPEN(70,FILE=fname,action='write',status='unknown',form='formatted')

      do i_z=1,shape_map(3)
       do i_y=1,shape_map(2)
        do i_x=1,shape_map(1)

        Write(70,'(3(I6,2x),4(f14.8,2x))') i_x,i_y,i_z,(map_vort(j,i_x,i_y,i_z),j=1,3), &
     & map_eul(i_x,i_y,i_z)

        enddo
       enddo
      enddo
      close(70)

      end subroutine topoplot_3D

      end module m_topoplot
