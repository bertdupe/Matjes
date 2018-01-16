      subroutine check_motif(motif)
      use m_derived_types
      implicit none
      type (cell), intent(inout) :: motif
      integer :: natom,nmag,i
      logical :: origin

      origin=.False.
      nmag=0
      i=1
      natom=size(motif%mom)

      do while (i.lt.(natom+1))
       if (abs(motif%mom(i)).gt.1.0d-5) then
        nmag=nmag+1
        motif%i_mom(i)=.True.
        else
        motif%i_mom(i)=.False.
       endif
      i=i+1
      enddo

      i=1
      do while (i.le.natom)
       if ((abs(sum(motif%pos(i,1:3))).lt.1.0d-8).and.(motif%i_mom(i))) origin=.True.
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

      subroutine update_lattice(motif,N_site)
      use m_derived_types
      implicit none
      type (cell), intent(in) :: motif
      integer, intent(inout) :: N_site

! number of magnetic atom per unit cell
      N_site=count(motif%i_mom)*N_site

      end subroutine update_lattice
