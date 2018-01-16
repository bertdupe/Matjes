!!!!!!!!!!!
! very easy routine that picks up a spin randomly
!!!!
      module m_choose_spin

      interface choose_spin

       module procedure choose_spin_serial
#ifdef CPP_MPI
       module procedure choose_spin_parallel
#endif

      end interface

      contains

      subroutine choose_spin_serial(Ilat,state,i_stone,mu_s,shape_spin,spin)
      use mtprng
      implicit none
      integer, intent(inout) :: Ilat(:)
      type(mtprng_state), intent(inout) :: state
      logical, intent(in) :: i_stone
      integer, intent(in) :: shape_spin(:)
      real(kind=8), intent(in) :: spin(:,:,:,:,:)
      real(kind=8), intent(out) :: mu_s
! internal variables
      real(kind=8) :: Choice
      integer :: i

      Choice=0.0d0

      Ilat=1

!     Choose a random spin place
      do i=1,3
#ifdef CPP_MRG
       choice=mtprng_rand_real1(state)
#else
       CALL RANDOM_NUMBER(Choice)
#endif
       if (shape_spin(i+1).ne.1) Ilat(i) = 1 + NINT(Choice*Dble(shape_spin(i+1))-0.5d0)
!       if (Ilat(i).gt.shape_spin(i)) Ilat(i)=1
      enddo

      if (shape_spin(5).gt.1) then
#ifdef CPP_MRG
       choice=mtprng_rand_real1(state)
#else
       CALL RANDOM_NUMBER(Choice)
#endif
       Ilat(4) = 1 + NINT(Choice*Dble(shape_spin(5))-0.5d0)
!       if (Ilat(4).gt.count(motif%i_m)) Ilat(4)=1
      endif

       if (i_stone) then
#ifdef CPP_MRG
        choice=mtprng_rand_real1(state)
#else
        CALL RANDOM_NUMBER(Choice)
#endif
        mu_s=spin(7,Ilat(1),Ilat(2),Ilat(3),Ilat(4))+(2.0d0*Choice-1.0d0)
       endif

      end subroutine choose_spin_serial

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef CPP_MPI
      subroutine choose_spin_parallel(Ilat,state,i_separate,i_average,i_ghost,i_stone,n_world, &
    &   mu_s,irank_working,shape_spin,MPI_COMM,spin,start)
      use mtprng
      implicit none
      integer, intent(inout) :: Ilat(:)
      type(mtprng_state), intent(inout) :: state
      integer, intent(in) :: n_world,start(:),N(:)
      logical, intent(in) :: i_separate,i_average,i_ghost,i_stone
<<<<<<< HEAD
      integer, intent(in) :: irank_working,shape_spin(:),MPI_COMM,n_world,start(:)
=======
>>>>>>> origin/Bertrand
      real(kind=8), intent(in) :: spin(:,:,:,:,:)
      real(kind=8), intent(out) :: mu_s
! internal variables
      real(kind=8) :: Choice
<<<<<<< HEAD
      integer :: ierr,N(3)
=======
      integer, intent(in) :: irank_working,shape_spin(:),MPI_COMM
      integer :: ierr
>>>>>>> origin/Bertrand
      integer :: i

      include 'mpif.h'
      Choice=0.0d0
      N=0
      N(1)=start(1)+shape_spin(2)
      N(2)=start(2)+shape_spin(3)
      N(3)=start(3)+shape_spin(4)

      if ((.not.i_separate).and.(.not.i_average).and.(.not.i_ghost)) then

        if (irank_working.eq.0) then
         do i=1,n_world
#ifdef CPP_MRG
          Choice=mtprng_rand_real1(state)
#else
          CALL RANDOM_NUMBER(Choice)
#endif
          Ilat(i) = 1 + NINT(Choice*Dble(N(i))-0.5d0)+start(i)

         enddo

         if (shape_spin(5).gt.1) then
#ifdef CPP_MRG
          Choice=mtprng_rand_real1(state)
#else
          CALL RANDOM_NUMBER(Choice)
#endif
          Ilat(4) = 1 + NINT(Choice*Dble(shape_spin(5))-0.5d0)

         endif

        if (i_stone) then
#ifdef CPP_MRG
          Choice=mtprng_rand_real1(state)
#else
          CALL RANDOM_NUMBER(Choice)
#endif
          mu_s=spin(7,Ilat(1),Ilat(2),Ilat(3),Ilat(4))+(2.0d0*Choice-1.0d0)
        else
         mu_s=spin(7,Ilat(1),Ilat(2),Ilat(3),Ilat(4))
        endif
       endif

       call MPI_BCAST(Ilat,4,MPI_INTEGER,0,MPI_COMM,ierr)
       if (i_stone) call MPI_BCAST(mu_s,1,MPI_REAL,0,MPI_COMM,ierr)

      else

         do i=1,n_world
#ifdef CPP_MRG
          Choice=mtprng_rand_real1(state)
#else
          CALL RANDOM_NUMBER(Choice)
#endif
          Ilat(i) = 1 + NINT(Choice*Dble(N(i))-0.5d0)+start(i)

         enddo

         if (shape_spin(5).gt.1) then
#ifdef CPP_MRG
          Choice=mtprng_rand_real1(state)
#else
          CALL RANDOM_NUMBER(Choice)
#endif
          Ilat(4) = 1 + NINT(Choice*Dble(shape_spin(5))-0.5d0)
         endif

         if (i_stone) then
#ifdef CPP_MRG
          Choice=mtprng_rand_real1(state)
#else
          CALL RANDOM_NUMBER(Choice)
#endif
          mu_s=spin(7,Ilat(1),Ilat(2),Ilat(3),Ilat(4))+(2.0d0*Choice-1.0d0)
          else
          mu_s=spin(7,Ilat(1),Ilat(2),Ilat(3),Ilat(4))
         endif
      endif

      end subroutine choose_spin_parallel
#endif

      end module m_choose_spin
