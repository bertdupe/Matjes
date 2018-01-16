      integer function qorienplot(map_qorien,my_lattice)
      use m_constants, only : k_b,pi
      use m_lattice, only : spin
      use m_parameters, only : kt
      use m_derived_types
      implicit none
      real(kind=8), dimension(:,:,:,:), intent(in) :: map_qorien
      type(lattice), intent(in) :: my_lattice
! internal
      integer :: i_x,i_y,i_z,j,i,dim_lat(3)
      real(kind=8) :: translation(3),normq,net(3,3)
      character(len=30) :: fname,toto

      dim_lat=my_lattice%dim_lat
      net=my_lattice%areal
      qorienplot=0
      translation=(net(1,:)+net(2,:)+net(3,:))/2.0d0

      write(6,'(a)') 'plot the topological charge'

      write(fname,'(f8.4)')kT/k_B
      toto=trim(adjustl(fname))
      write(fname,'(a,18a,a)')'qorienmap',(toto(i:i),i=1, &
       len_trim(toto)),'.dat'
      OPEN(70,FILE=fname,action='write',status='unknown',form='formatted')
       do i_z=1,dim_lat(3)
        do i_y=1,dim_lat(2)
         do i_x=1,dim_lat(1)

        normq=sqrt(map_qorien(1,i_x,i_y,i_z)**2+map_qorien(2,i_x,i_y,i_z)**2+map_qorien(3,i_x,i_y,i_z)**2)
        if (normq.lt.1.0d-7) normq=1.0d0
        Write(70,'(8(E20.10E3,2x))') (Spin(j,i_x,i_y,i_z,1)+translation(j),j=1,3), &
     &   acos(map_qorien(3,i_x,i_y,i_z)/normq),atan2(map_qorien(2,i_x,i_y,i_z),map_qorien(1,i_x,i_y,i_z)), &
     &   (map_qorien(j,i_x,i_y,i_z)/pi(4.0d0),j=1,3)

         enddo
        enddo
       enddo
      close(70)

      end function qorienplot

