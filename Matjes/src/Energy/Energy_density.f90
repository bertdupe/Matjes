! =====================================================================
! Subroutine that gives the energy density on a lattice, given an input
! file containing the lattice.
! Author : Charles Paillard
! Starting Writing date : 21.05.2013 (May 21st, 2013)
! Last modification : May, 28th 2013
! Version : 1.1
! =====================================================================

      SUBROUTINE EnergyDensity(EnOutput,EnDensityOutput,spin,shape_spin,tableNN,shape_tableNN,masque,shape_masque,indexNN,shape_index,h_ext)
      use m_rw_lattice, only : net,world,motif
      use m_parameters, only : J_ij,EA,J_il,i_biq,i_four,i_DM,i_dip,DM,Periodic_log,J_z
      use m_constants, only : pi
      use m_energy
      use m_table_dist
      Implicit None
! input
      integer, intent(in) :: shape_index(2),shape_spin(5),shape_tableNN(6),shape_masque(4)
      real(kind=8), intent(in) :: spin(shape_spin(1),shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5)),h_ext(3)
      integer, intent(in) :: tableNN(shape_tableNN(1),shape_tableNN(2),shape_tableNN(3),shape_tableNN(4),shape_tableNN(5),shape_tableNN(6))
      integer, intent(in) :: masque(shape_masque(1),shape_masque(2),shape_masque(3),shape_masque(4))
      integer, intent(in) :: indexNN(shape_index(1),shape_index(2))
      Character (LEN = *), intent(inout) ::EnOutput, EnDensityOutput
! internals
      real(kind=8), allocatable ::EnDenTab(:)
      Integer::ind,j,nb_col,NColEnergy,NStartDensity,col,i_x,i_y,i_z,i_m
      real (kind=8), allocatable :: tabledist(:,:)
      logical::dexists
      real(kind=8) :: SNeigh, Sprevious
      integer :: N_Nneigh,Nei_z,N_nei_exch,N_nei_dm,Nei_DM,Nei_il
      integer :: alloc_check,nmag
! Table containing all energy terms. The different columns contain the
! following information :
! --> 1  = X position on the lattice
! --> 2  = Y position on the lattice
! --> 3  = Exchange energy
! --> 4  = Dipole-Dipole Interaction energy
! --> 5  = Zeeman energy
! --> 6  = Dzyaloshinskii-Moriya interaction energy
! --> 7  = Anisotropy energy
! --> 8  = Biquadratic interaction energy
! --> 9  = Four spin interaction energy
! --> 10 = Sum of all these contributions, that we will call total
!          energy
! --> 11-(10+N_Nneigh) : Total energy to the point due to the k-th
! neighbour
! --> (11+N_Nneigh) - (10+2*N_Nneigh) : surface energy density due to the
! k-th neighbours

! initialization of the local variables
      N_nei_exch=2
      N_nei_dm=2
!! locally change the motif for density of energy

      N_Nneigh=max(size(abs(J_ij),1),size(J_il))
      Nei_il=size(J_il)
      Nei_z=size(J_z)
      Nei_DM=count(abs(DM).gt.1.0d-8)
      nmag=shape_spin(5)

! take into account the number of neighbors in case of none periodic boundary condition
      if (.not.all(Periodic_log)) then
       N_nei_exch=sum(indexNN(1:N_Nneigh,1))
       N_nei_dm=sum(indexNN(1:Nei_DM,1))
      endif

      if ((nmag.eq.1).and.(shape_spin(5).eq.1)) then
       allocate(tabledist(N_Nneigh,1),stat=alloc_check)
       if (alloc_check.ne.0) stop 'can not allocate table of distance'
       tabledist=0.0d0
       call Tdist(net,N_Nneigh,world,1,motif,tabledist(:,1))
      else
       allocate(tabledist(N_Nneigh,2),stat=alloc_check)
       if (alloc_check.ne.0) stop 'can not allocate table of distance'
       tabledist=0.0d0
       call Tdist(net,N_Nneigh,Nei_z,Nei_il,world,2,motif,tabledist(:,:))
      endif

      nb_col   = 10 + 2*N_Nneigh+7
      EnOutput = TRIM(ADJUSTL(EnOutput))
      EnDensityOutput = TRIM(ADJUSTL(EnDensityOutput))
      Sneigh   = 1.0d0
      Sprevious= 1.0d0

      Open(unit = 666,file=EnOutput,form='formatted',action ='write')
      Open(unit = 667,file=EnDensityOutput,form='formatted',action='write')

      Allocate(EnDenTab(nb_col),stat=alloc_check)
      if (alloc_check.ne.0) stop 'can not allocate energy density matrix'
! If you change the number of Near neighbours, you might have some
! troubles with the format 1998. Please, always write it in the
! following form : format( (10+2*N_Nneigh)f14.8 ).
99999  format(41f14.8)
99998  format(33f14.8)
99997  format(21f14.8)


      do i_z=1,shape_spin(4)
       do i_y=1,shape_spin(3)
        do i_x=1,shape_spin(2)
        EnDenTab=0.0d0
        EnDenTab(1)=Spin(1,i_x,i_y,i_z,1)
        EnDenTab(2)=Spin(2,i_x,i_y,i_z,1)

        do i_m=1,shape_spin(5)
!           do i_m=1,1

      EnDenTab(3)=EnDenTab(3)+Exchange(i_x,i_y,i_z,i_m,spin,shape_spin,tableNN,masque,indexNN)
!      EnDenTab(3)=EnDenTab(3)+Exchange(i_x,i_y,i_z,i_m)/dble(count(masque(2:N_nei_exch+1,i_x,i_y,i_z).ne.0))
#ifdef CPP_BRUTDIP
      if (i_dip) EnDenTab(4)=EnDenTab(4)+dipole(i_x,i_y,i_z,i_m,spin,shape_spin,masque)
#else
      if (i_dip) EnDenTab(4)=EnDenTab(4)+fftdip(i_x,i_y,i_z,i_m)
#endif
      EnDenTab(5)=EnDenTab(5)+Zeeman(i_x,i_y,i_z,i_m,spin,shape_spin,masque,h_ext)
      if (i_DM) EnDenTab(6)=EnDenTab(6)+DMenergy(i_x,i_y,i_z,i_m,spin,shape_spin,tableNN,masque,indexNN)

!      if (i_DM) EnDenTab(6)=EnDenTab(6)+DMenergy(i_x,i_y,i_z,i_m)/dble(count(masque(2:N_nei_dm+1,i_x,i_y,i_z).ne.0))
      EnDenTab(7)=EnDenTab(7)+anisotropy(i_x,i_y,i_z,i_m,EA,spin,shape_spin,masque)
      if (i_biq) EnDenTab(8)=EnDenTab(8)+biquadratic(i_x,i_y,i_z,i_m,spin,shape_spin,tableNN,masque,indexNN)
      if (i_four) EnDenTab(9)=EnDenTab(9)+fourspin(i_x,i_y,i_z,i_m,spin,shape_spin,masque)
        enddo

      EnDenTab(10)= sum(EnDenTab(3:9))

       Do j=1, N_Nneigh

            SNeigh=pi(1.0d0)*tabledist(j,1)**2

! Filling the Zeeman energy density, the DM energy density, The
! anisotropy energy density, the biquadratic energy density, the
! fourspin energy density and the total energy density
            EnDenTab(10+2*N_Nneigh+1) = EnDenTab(5)/SNeigh
            EnDenTab(10+2*N_Nneigh+2) = EnDenTab(6)/SNeigh
            EnDenTab(10+2*N_Nneigh+3) = EnDenTab(7)/SNeigh
            EnDenTab(10+2*N_Nneigh+4) = EnDenTab(8)/SNeigh
!The density of four spin energy has not been implemented yet.
            EnDenTab(10+2*N_Nneigh+5) = EnDenTab(4)/SNeigh
            EnDenTab(10+2*N_Nneigh+6) = EnDenTab(9)/SNeigh
            EnDenTab(10+2*N_Nneigh+7) = EnDenTab(10)/SNeigh

            if (j.eq.1) then
            
            EnDenTab(10+j) = ExchN(i_x,i_y,i_z,1, j)
            EnDenTab(10+N_Nneigh+j) = EnDenTab(10+j)/SNeigh
            Sprevious = SNeigh
!      WRITE(*,*) 'Energy_density l.71 OK boucle numero ', ind, ' , ', j
!       STOP

            else

            EnDenTab(10+j) = ExchN(i_x,i_y,i_z,1, j)-ExchN(i_x,i_y,i_z,1, j-1)
            EnDenTab(10+N_Nneigh+j) = EnDenTab(10+j)&
                                  /(SNeigh-Sprevious)
            Sprevious = SNeigh
!      WRITE(*,*) 'Energy_density l.71 OK boucle numero ', ind, ' , ', j
!      STOP

            end if

         EndDo

         
      NColEnergy   = 10+N_Nneigh
      NStartDensity = 10+N_Nneigh+1
      Write(666,99998) (EnDenTab(col), col=1,NColEnergy)
      Write(667,99997) (EnDenTab(col), col=1,2), (EnDenTab(col), col=  &
                         NStartDensity, nb_col)
!      Write(666,99999) (EnDenTab(col), col=1,40)

        enddo
       enddo
      enddo

      Deallocate(EnDenTab,tabledist)

      Close(666)
      Close(667)

      END SUBROUTINE EnergyDensity


