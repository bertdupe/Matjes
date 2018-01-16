! ===============================================================
      SUBROUTINE InitSpin(my_lattice,motif)
      use m_get_position
      use m_init_modes
      use m_init_config
      use m_vector, only : cross,norm
      use m_constants
      use m_derived_types
      use m_get_random
      Implicit none
! variables that come in
      type (cell), intent(in) :: motif
      type(lattice), intent(inout) :: my_lattice
!     Slope Indexes for three dim spins
      LOGICAL :: i_exi,i_init_config
!     Absolute value of a spin
      integer, parameter  :: io=9
      integer :: dim_lat(3),nmag
      integer :: i_x,i_y,i_z,i_m
      real(kind=8) :: r(3,3),mu_S

#ifdef CPP_MPI
      include 'mpif.h'
      if (irank.eq.0) then 
#endif
      mu_S=motif%mom(1)
      dim_lat=my_lattice%dim_lat
      r=my_lattice%areal

!     Check if spin lattice input is present
      inquire(file='init-modes.dat',exist=i_exi)

! part that reads the old lattice
      if (i_exi) then
       write(6,'(a)') 'Read spin lattice from SpinSTMi.dat'
       call get_init_modes('spinSTMi.dat',my_lattice,motif)
       return

      endif

! initialize the lattice

!truc bizare
!        do i=1,dim_lat(1)
!        do j=1,dim_lat(2)
!         Spin(4:6,i,j,1,1)=Spin(4:6,i,72,1,1)
!         Spin(4:6,i,j,1,2)=Spin(4:6,i,72,1,2)
!
!        enddo
!        enddo
! Check for reordering the lattice

      i_init_config=.False.
      inquire(file='init.config',exist=i_init_config)
      if (i_init_config) then
         call init_config('init.config',my_lattice,motif)
         return
      endif

      nmag=count(motif%i_mom)
      do i_m=1,nmag
         do i_z=1,dim_lat(3)
            do i_y=1,dim_lat(2)
               do i_x=1,dim_lat(1)

! fix the spin direction as random

                my_lattice%l_modes(i_x,i_y,i_z,i_m)%w(:)=get_rand_classic(3,1.0d0)

               enddo
            enddo
         enddo
      enddo


#ifdef CPP_MPI
      endif

      spin=bcast(Spin,7,dim_lat(1),dim_lat(2),dim_lat(3),count(motif%i_m),MPI_REAL8,0,MPI_COMM)

#endif

#ifdef CPP_DEBUG
      do i_m=1,nmag
       do i_z=1,dim_lat(3)
        do i_y=1,dim_lat(2)
         do i_x=1,dim_lat(1)
          write(*,*) my_lattice%l_modes(i_x,i_y,i_z,i_m)%w(:)
         enddo
        enddo
       enddo
      enddo
      stop
#endif

      END SUBROUTINE InitSpin
! ===============================================================


