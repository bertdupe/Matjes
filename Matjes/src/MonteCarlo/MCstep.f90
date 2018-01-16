!
! ===============================================================
!
      SUBROUTINE MCstep(state,E_total,E,Magnetization,kt,acc,rate,nb,cone, &
    &  spin,shape_spin,tableNN,shape_tableNN,masque,shape_masque,indexNN,shape_index,h_ext,EA, &
    & i_biq,i_dip,i_DM,i_four,i_stone,ising,equi,overrel,sphere,underrel,my_lattice)
      use m_derived_types
      use m_sampling
      use m_choose_spin
      use m_relaxtyp
      use m_createspinfile
      use m_local_energy
#ifdef CPP_MPI
      use m_parameters, only : d_mu,n_ghost,i_separate,i_average,i_ghost
      use m_make_box, only : ghost_border,rank_nei
      use m_mpi_prop, only : isize,MPI_COMM,N,start,MPI_COMM_BOX,irank_working,irank_box
#else
      use m_parameters, only : d_mu
#endif
      use m_constants
      use mtprng
#ifndef CPP_BRUTDIP
      use m_setup_dipole, only : mmatrix
#endif
      Implicit none
! input
      type(lattice), intent(in) :: my_lattice
      integer, intent(in) :: shape_index(2),shape_spin(5),shape_tableNN(6),shape_masque(4)
      integer, intent(in) :: tableNN(shape_tableNN(1),shape_tableNN(2),shape_tableNN(3),shape_tableNN(4),shape_tableNN(5),shape_tableNN(6))
      integer, intent(in) :: masque(shape_masque(1),shape_masque(2),shape_masque(3),shape_masque(4))
      integer, intent(in) :: indexNN(shape_index(1),shape_index(2))
      logical, intent(in) :: i_biq,i_dip,i_DM,i_four,i_stone,ising,equi,overrel,sphere,underrel
      real(kind=8), intent(in) :: kt,h_ext(3),EA(3)
      real(kind=8), intent(inout) :: spin(shape_spin(1),shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5))
      type(mtprng_state), intent(inout) :: state
      real(kind=8), intent(inout) :: E_total,Magnetization(3),E(8),acc,rate,cone,nb
!     Random spin coordinate
      Integer :: Ilat(4)
! magnetic moments
      real(kind=8) :: mu_s,mu_old
!     Energy difference Delta E and -DE/kT, Dmag and Dq
      real(kind=8) :: DE,E_new,E_old,Dmag(3)
      real(kind=8) :: tmp
!     memory value for Theta
      real(kind=8) :: choice
!     flipped Spin
      real(kind=8) :: S_new(3),S_old(3)
      logical :: accept
! decomposition of the energy
      real(kind=8) :: E_dec_new(8),E_dec_old(8)


#ifdef CPP_MPI
      real(kind=8) :: mpi_DE(isize),mpi_S_store(3,isize),mpi_mu(isize)
      integer :: imin(1),irank_discuss,irank_send
! shall the comm occur
      logical :: discuss,at_border,discuss_trans

! control of the point to point comm
      integer :: reqs(4)

      include 'mpif.h'

      mpi_DE=0.0d0
      mpi_S_store=0.0d0
      mpi_mu=0.0d0
      discuss=.False.
      at_border=.False.
      discuss_trans=.False.
      irank_discuss=0
      irank_send=0

#endif

      Ilat=1
      E_dec_new=0.0d0
      E_dec_old=0.0d0

#ifdef CPP_MPI
      call choose_spin(Ilat,state,i_separate,i_average,i_ghost,i_stone,n_world, &
    &   mu_s,irank_working,shape_spin,MPI_COMM,spin,start)
#else
      call choose_spin(Ilat,state,i_stone,mu_s,shape_spin,spin)
#endif

! check if the spin exists
      if (masque(1,Ilat(1),Ilat(2),Ilat(3)).eq.0) return
!     Initialize Store variables
      S_new   =0.0d0

! cone angle update and maximal magnetic moment change
        if ((rate.gt.0.5d0).and.(cone.lt.pi(1.0d0)))then
           cone=cone+0.0001d0
        elseif ((rate.lt.0.50d0).and.(cone.gt.0.01d0)) then
           cone=cone-0.0001d0
        endif
        if (i_stone) then
         if ((rate.gt.0.5d0).and.(d_mu.lt.1.0d0))then
           d_mu=d_mu+0.0001d0
         elseif ((rate.lt.0.50d0).and.(d_mu.gt.1.0d-5)) then
           d_mu=d_mu-0.0001d0
         endif
        endif
!---------------------------------------
! here are the different sampling
! first the sphere sampling
!---------------------------------------

      if (sphere) then
        S_new=sphereft(Spin(4:6,Ilat(1),Ilat(2),Ilat(3),Ilat(4)),cone,state)
!---------------------------------------
! Algorithm square sampling
!---------------------------------------
      elseif (equi) then
        S_new=equirep(Spin(4:6,Ilat(1),Ilat(2),Ilat(3),Ilat(4)),cone,state)
      endif
!---------------------------------
! different relaxation process
!---------------------------------
      if (underrel) then
        S_new=underrelax(Spin(4:6,Ilat(1),Ilat(2),Ilat(3),Ilat(4)),Ilat(1),Ilat(2),Ilat(3),Ilat(4),spin,shape_spin,indexNN, &
              &    shape_index,masque,shape_masque,tableNN,shape_tableNN,h_ext,my_lattice)
      elseif (overrel) then
        S_new=overrelax(Spin(4:6,Ilat(1),Ilat(2),Ilat(3),Ilat(4)),Ilat(1),Ilat(2),Ilat(3),Ilat(4),spin,shape_spin,indexNN, &
              &    shape_index,masque,shape_masque,tableNN,shape_tableNN,h_ext,my_lattice)
      endif

      if (ising) S_new=-Spin(4:6,Ilat(1),Ilat(2),Ilat(3),Ilat(4))
!----------------------------------
!       Calculate the energy difference if this was flipped
!       and decider, if the Spin flip will be performed
!----------------------------------
!Energy of old configuration
        E_dec_old=local_energy_MC(i_DM,i_four,i_biq,i_dip,i_stone,EA,Ilat, &
          & spin,shape_spin,tableNN,shape_tableNN,masque,shape_masque,indexNN,shape_index,h_ext,my_lattice)

        E_old=sum(E_dec_old)

!        write(*,*) 'energy'
!        write(*,*) Exchange(Ilat(1),Ilat(2),Ilat(3),Ilat(4))
!        write(*,*) Zeeman(Ilat(1),Ilat(2),Ilat(3),Ilat(4))
!        if (i_DM) write(*,*) DMenergy(Ilat(1),Ilat(2),Ilat(3),Ilat(4))
!        write(*,*) anisotropy(Ilat(1),Ilat(2),Ilat(3),Ilat(4),EA)
!        if (i_dip) write(*,*) dipole(Ilat)
!        pause
!Energy of the new configuration
        S_old=Spin(4:6,Ilat(1),Ilat(2),Ilat(3),Ilat(4))
        if (i_stone) mu_old=Spin(7,Ilat(1),Ilat(2),Ilat(3),Ilat(4))

        Spin(4:6,Ilat(1),Ilat(2),Ilat(3),Ilat(4))=S_new
        if (i_stone) Spin(7,Ilat(1),Ilat(2),Ilat(3),Ilat(4))=mu_s

        E_dec_new=local_energy_MC(i_DM,i_four,i_biq,i_dip,i_stone,EA,Ilat, &
          & spin,shape_spin,tableNN,shape_tableNN,masque,shape_masque,indexNN,shape_index,h_ext,my_lattice)

        E_new=sum(E_dec_new)

        Spin(4:6,Ilat(1),Ilat(2),Ilat(3),Ilat(4))=S_old
        if (i_stone) Spin(7,Ilat(1),Ilat(2),Ilat(3),Ilat(4))=mu_old

#ifndef CPP_BRUTDIP
         if (i_dip) mmatrix(:,Ilat(1),Ilat(2),Ilat(3))=S_old*Spin(7,Ilat(1),Ilat(2),Ilat(3),Ilat(4))
#endif

!! variation of the energy for this step
        DE=E_new-E_old

!---------------------------------------------------
! Here is the different way to compute de next 
! candidate for the spin
!---------------------------------------------------

! End of the boradcast
#ifdef CPP_MPI
      if ((.not.i_separate).and.(.not.i_average).and.(.not.i_ghost)) then

        call mpi_gather(DE,1,MPI_REAL8,mpi_DE,1,MPI_REAL8,0,MPI_COMM,ierr)
        call mpi_gather(S_new,3,MPI_REAL8,mpi_S_store(:,:),3,MPI_REAL8,0,MPI_COMM,ierr)
        if (i_stone) call mpi_gather(mu_s,1,MPI_REAL8,mpi_mu,3,MPI_REAL8,0,MPI_COMM,ierr)

         if (irank_working.eq.0) then
          imin=minloc(mpi_DE)
          DE=mpi_DE(imin(1))

!           do i_store=1,3
!           S_store(i_store)=mpi_S_store(3*(imin(1)-1)+i_store)
!           enddo
          S_new=mpi_S_store(:,imin(1))
          if (i_stone) mu_s=mpi_mu(imin(1))
         endif

        call mpi_bcast(DE,1,MPI_REAL8,0,MPI_COMM,ierr)
        call mpi_bcast(S_new,3,MPI_REAL8,0,MPI_COMM,ierr)
        if (i_stone) call mpi_bcast(mu_s,1,MPI_REAL8,0,MPI_COMM,ierr)

       endif
#endif

        nb=nb+1.0d0
        accept=.False.
! security in case kt is 0
        if (kt.gt.1.0d-10) then
!       Calculate the probability of this flip
         tmp=-DE/kT
         if (DE.lt.0.0d0) tmp=2.0d0
        else
         tmp=-100.0d0
!       Accept the change of the new energy is lower than the old one
         if (DE.lt.0.0d0) tmp=2.0d0
        endif

        if (exp(tmp).gt.1.0d0) then
         accept=.True.
         else
#ifdef CPP_MPI
#ifdef CPP_MRG
          choice=mtprng_rand_real1(state)
#else
          CALL RANDOM_NUMBER(Choice)
#endif
          if ((.not.i_separate).and.(.not.i_average).and.(.not.i_ghost)) call mpi_bcast(Choice,1,MPI_REAL8,0,MPI_COMM,ierr)
#else
#ifdef CPP_MRG
          choice=mtprng_rand_real1(state)
#else
          CALL RANDOM_NUMBER(Choice)
#endif
#endif
          if (Choice.lt.exp(tmp)) Then
           accept=.True.
          endif
        endif

        if (accept) then

        Dmag=-Spin(4:6,Ilat(1),Ilat(2),Ilat(3),Ilat(4))+S_new

! update the spin
        Spin(4:6,Ilat(1),Ilat(2),Ilat(3),Ilat(4))=S_new

! update the quantities
        E_total=E_total+DE
        E=E+E_dec_new-E_dec_old
        Magnetization=Magnetization+Dmag

        if (i_stone) Spin(7,Ilat(1),Ilat(2),Ilat(3),Ilat(4))=mu_s
#ifndef CPP_BRUTDIP
        if (i_dip) mmatrix(:,Ilat(1),Ilat(2),Ilat(3))=S_new*Spin(7,Ilat(1),Ilat(2),Ilat(3),Ilat(4))
#endif
       acc=acc+1.0d0
       endif

#ifdef CPP_MPI
! if the move is accepted, one has to make sure that the updated spin does not lie at one interface

        if (i_ghost) then
         at_border=ghost_border(1,Ilat(1),Ilat(2),Ilat(3))
!         write(*,*) at_border ,Ilat(1),Ilat(2), irank_working, irank_box
         do i=1,n_ghost
          if ((accept).and.(at_border).and.((i-1).eq.irank_box)) then
           discuss_trans=.True.
           irank_send=irank_box
          endif
         call mpi_allreduce(discuss_trans,discuss,1,mpi_logical,MPI_LOR,MPI_COMM_BOX,ierr)
         call mpi_allreduce(irank_send,irank_discuss,1,mpi_integer,MPI_sum,MPI_COMM_BOX,ierr)

! if it lies at one of the interface then update the spins of the neighbors

         if (discuss) call MC_transfer(Ilat,S_new,irank_discuss,spin,shape_spin)

         discuss_trans=.False.
         irank_send=0
         enddo
        endif
#endif

      rate=acc/nb
!      write(*,*) acc,nb,kt/k_B,irank_working
! part that does the parallel tempering

      END SUBROUTINE MCstep
! ===============================================================
