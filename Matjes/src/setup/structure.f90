      subroutine check_motif(motif)
      use m_derived_types
      implicit none
      type (cell), intent(inout) :: motif
      integer :: natom,nmag
      integer :: i,j
      logical :: origin

      origin=.False.
      nmag=0
      i=1
      natom=size(motif%mom)

      do while (i.lt.(natom+1))
       if (abs(motif%mom(i)).gt.1.0d-5) then
        nmag=nmag+1
        motif%i_m(i)=.True.
        else
        motif%i_m(i)=.False.
       endif
      i=i+1
      enddo

      i=1
      do while (i.le.natom)
       do j=1,size(motif%i_m)
        if ((motif%i_m(i)).and.(dabs(motif%mom(j)-motif%mom(i)).lt.1.0d-8)) then
         motif%type(i)=1
        else
         motif%type(i)=0
        endif
       enddo
      i=i+1
      enddo

      i=1
      do while (i.le.natom)
       if ((abs(sum(motif%pos(i,1:3))).lt.1.0d-8).and.(motif%i_m(i))) origin=.True.
      i=i+1
      enddo

      if (.not.origin) then
       write(6,*) "a magnetic atom must be at (0,0,0) or the neighbors indexation will not work"
       write(6,*) "I will stop"
       stop
      endif

#ifdef CPP_DEBUG
      i=1
      do while (i.lt.(natom+1))
       write(6,'(4f8.4)') (motif%pos(i,j),j=1,3),motif%mom(i)
       i=i+1
      enddo 
#endif

      end subroutine check_motif

      subroutine check_net(net,dim_lat)
      use m_vector, only : TripleProduct,norm
      implicit none
! net must be direct or the topo charge will be opposite... or not. I do not really know actually
! inout variable
      real(kind=8), intent(inout) :: net(3,3)
      integer, intent(in) :: dim_lat(3)
! dummy
      real(kind=8) :: dumnet(3,3)

      if (norm(net(1,:)).lt.0.00001d0) then
       write(6,'(a)') 'the x axis must be different from 0'
       write(6,'(a)') 'rotate your cell'
       stop
      endif

      if (TripleProduct(net(1,:),net(2,:),net(3,:)).lt.0.0d0) then
       write(*,*) 'indirect lattice found, several things are not compatible with that'
       write(*,*) 'please unter a lattice so that u.(vxw)>0'
       stop
      endif

      end subroutine check_net

      subroutine update_lattice(motif,N_site)
      use m_derived_types
      implicit none
      type (cell), intent(in) :: motif
      integer, intent(inout) :: N_site
      integer :: i,j

! number of magnetic atom per unit cell
      N_site=count(motif%i_m)*N_site

      end subroutine update_lattice
