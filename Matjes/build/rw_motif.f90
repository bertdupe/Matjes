      subroutine rw_motif(my_motif,my_lattice)
      use m_derived_types
      implicit none
! intent(inout)
      type(lattice), intent(inout) :: my_lattice
      type(cell), intent(out) :: my_motif
! internal
      integer, parameter :: io=2
      character(len=100) :: str
! check the allocation of memory
      integer :: alloc_check
      character(len=10) :: dummy
      integer :: fin,i,j,N_site,natom,dimension(3)

      dimension=my_lattice%dim_lat

      open (io,file='input',form='formatted',status='old',action='read')

      rewind(io)
      do
      read (io,'(a)',iostat=fin) str
       if (fin /= 0) exit
       str= trim(adjustl(str))
       if (len_trim(str)==0) cycle

      if (str(1:1) == '#' ) cycle
! We start to read the input

       if ( str(1:5) == 'motif') then
        backspace(io)
         read(io,*) dummy, natom, dummy
         allocate(my_motif%pos(natom,3),my_motif%mom(natom),my_motif%i_mom(natom),stat=alloc_check)
         if (alloc_check.ne.0) write(6,'(a)') 'can not allocate motif'
          my_motif%pos=0.0d0
          my_motif%mom=0.0d0
          my_motif%i_mom=.False.
          do i=1,natom
           read(io,*) (my_motif%pos(i,j),j=1,3),my_motif%mom(i)
           if (my_motif%mom(i).ne.0.0d0) my_motif%i_mom(i)=.True.
          enddo
       endif

      enddo
      close(io)


! size of the world
        if ((dimension(3).eq.1).and.(dimension(2).eq.1)) then
         allocate(my_lattice%world(1))
         my_lattice%world(1)=dimension(1)
         my_lattice%n_system=1
          if (count(my_motif%i_mom).gt.1) my_lattice%n_system=12
        elseif (dimension(3).eq.1) then
         allocate(my_lattice%world(2))
         my_lattice%world=(/dimension(1),dimension(2)/)
         my_lattice%n_system=2
          if (count(my_motif%i_mom).gt.1) my_lattice%n_system=22
        elseif ((dimension(3).eq.1).and.(dimension(2).eq.1).and.(dimension(1).eq.1)) then
         write(6,*) "dimension of the problem not correct"
         stop
        else
         allocate(my_lattice%world(3))
         my_lattice%world=(/dimension(1),dimension(2),dimension(3)/)
         my_lattice%n_system=3
         if (count(my_motif%i_mom).gt.1) my_lattice%n_system=32
        endif

      call check_motif(my_motif)

! takes into account that we might have more than one magnetic atom per unit cell

      call update_lattice(my_motif,N_site)

      end subroutine rw_motif
