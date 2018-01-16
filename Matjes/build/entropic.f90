      subroutine entropic(N_cell, &
            &    spin,shape_spin,tableNN,shape_tableNN,masque,shape_masque,indexNN,shape_index, &
            &    i_biq,i_dip,i_DM,i_four,i_stone,ising,equi,sphere,overrel,underrel, &
            &    EA,h_ext, &
            &    kTfin,kTini,n_Tsteps,my_lattice)
      use mtprng
      use m_choose_spin
      use m_local_energy
      use m_constants, only : pi,k_B
      use m_relaxtyp
      use m_sampling
      use m_average_MC
      use m_derived_types
      implicit none
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! variable part
      integer, intent(in) :: N_cell,n_Tsteps
      integer, intent(in) :: shape_spin(5),shape_tableNN(6),shape_masque(4),shape_index(2)
      integer, intent(in) :: tableNN(shape_tableNN(1),shape_tableNN(2),shape_tableNN(3),shape_tableNN(4),shape_tableNN(5),shape_tableNN(6))
      integer, intent(in) :: masque(shape_masque(1),shape_masque(2),shape_masque(3),shape_masque(4))
      integer, intent(in) :: indexNN(shape_index(1),shape_index(2))
      logical, intent(in) :: i_biq,i_dip,i_DM,i_four,i_stone,ising,equi,sphere,overrel,underrel
      real(kind=8), intent(in) :: EA(3),h_ext(3),kTfin,kTini
      type(lattice), intent(in) :: my_lattice

! inout
      real(kind=8), intent(inout) :: spin(shape_spin(1),shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5))
!------------------------------------------
! internal variable
      type(mtprng_state) :: state
      integer :: Ilat(4),n_steps,i_step,n_sweep,i,i_loop,n_loop
      real(kind=8) :: mu_s,cone,choice
      logical :: accept,i_exit
! coordinate of the spin in internal coordinate
      integer :: i_x_s,i_y_s,i_z_s,i_m_s
      real(kind=8) :: S_new(3),S_old(3)
! energy variable
      real(kind=8) :: E_min,E_max,E_delta,E_test
      real(kind=8) :: E_new,E_old,E_total,E_decompose(8)
      integer,allocatable :: count_E(:)
! histogram variables
      integer :: n_histo,n_old,n_new,min_histo,i_histo
      real(kind=8),allocatable :: lng(:)
      real(kind=8) :: lnf,lng_new,lng_old,delta_lng,flatness,mean_histo
! thermodynamical quantitites
      real(kind=8) :: energy(n_Tsteps),Z(n_Tsteps),temperatures(n_Tsteps),free_E(n_Tsteps), &
     &   entropy(n_Tsteps),energy_sq(n_Tsteps),heat_cap(n_Tsteps)
      real(kind=8) :: kT,DProb,DE
      integer :: i_T
#ifdef CPP_DEBUG
      integer :: vector(N_cell+1),length
      vector=0
      length=int(sum(Spin(6,:,:,:,:)))
#endif
! initialize variables
      Ilat=1
      mu_s=0.0d0
      E_min=0.0d0
      E_max=0.0d0
      E_new=0.0d0
      E_old=0.0d0
      E_delta=0.0d0
      n_histo=82
      allocate(lng(n_histo))
      energy=0.0d0
      heat_cap=0.0d0
      energy_sq=0.0d0
      entropy=0.0d0
      free_E=0.0d0
      Z=0.0d0
      temperatures=0.0d0
      lng=0.0d0
      cone=pi(1.0d0)
      S_new=0.0d0
      S_old=0.0d0
      n_steps=100000000
      n_loop=20
      lnf=0.50d0
      lng_new=0.0d0
      lng_old=0.0d0
      delta_lng=0.0d0
      n_old=0
      n_new=0
      accept=.False.
      allocate(count_E(n_histo))
      count_E=0
      E_test=0.0d0
      flatness=0.9d0
      i_exit=.False.

! initializing the variables above
      call DeriveValue(N_cell,spin,shape_spin,tableNN,shape_tableNN,masque,shape_masque,indexNN,shape_index,E_decompose, &
           &   i_four,i_dm,i_biq,i_dip,i_stone,EA,h_ext)
      E_total=sum(E_decompose)

      write(6,'(a)') "in entropic sampling subroutine"
! find the energy range

      E_min=-8.0d-3*real(N_cell)
      E_max=8.0000001d-3*real(N_cell)
#ifdef CPP_DEBUG
      if ((E_total-E_min).lt.0.0d0) then
        write(6,'(a,f12.7)') 'the minimum of energy ', E_min
        write(6,'(a,f12.7)') 'is smaller than the energy ', E_total
        stop
      endif
#endif

      E_delta=(E_max-E_min)/real(n_histo)
      write(6,'(a,E20.10E3)') 'The size of the energy bin is ', E_delta

! open the results file

      OPEN(7,FILE='histogram.dat',action='write',status='unknown',position='append',form='formatted')

      do i_loop=1,n_loop
          do i_step=1,n_steps
! pic a spin
             call choose_spin(Ilat,i_stone,mu_s,shape_spin,spin)

             i_x_s=Ilat(1)
             i_y_s=Ilat(2)
             i_z_s=Ilat(3)
             i_m_s=Ilat(4)

             S_old=Spin(4:6,i_x_s,i_y_s,i_z_s,i_m_s)
!---------------------------------------
! here are the different sampling
! first the sphere sampling
!---------------------------------------
             if (sphere) then
                S_new=sphereft(Spin(4:6,i_x_s,i_y_s,i_z_s,i_m_s),cone)
!---------------------------------------
! Algorithm square sampling
!---------------------------------------
             elseif (equi) then
                S_new=equirep(Spin(4:6,i_x_s,i_y_s,i_z_s,i_m_s),cone)
             endif
!---------------------------------
! different relaxation process
!---------------------------------
             if (underrel) then
                S_new=underrelax(Spin(4:6,i_x_s,i_y_s,i_z_s,i_m_s),i_x_s,i_y_s,i_z_s,i_m_s,spin,shape_spin,indexNN, &
                      &  shape_index,masque,shape_masque,tableNN,shape_tableNN,h_ext,my_lattice)
             elseif (overrel) then
                S_new=overrelax(Spin(4:6,i_x_s,i_y_s,i_z_s,i_m_s),i_x_s,i_y_s,i_z_s,i_m_s,spin,shape_spin,indexNN, &
                      &  shape_index,masque,shape_masque,tableNN,shape_tableNN,h_ext,my_lattice)
             endif

             if (ising) S_new=-Spin(4:6,i_x_s,i_y_s,i_z_s,i_m_s)

! calculate the energy before and after

             E_decompose=local_energy_MC(i_DM,i_four,i_biq,i_dip,i_stone,EA,Ilat, &
      & spin,shape_spin,tableNN,shape_tableNN,masque,shape_masque,indexNN,shape_index,h_ext,my_lattice)
             E_old=sum(E_decompose)

             Spin(4:6,i_x_s,i_y_s,i_z_s,i_m_s)=S_new

             E_decompose=local_energy_MC(i_DM,i_four,i_biq,i_dip,i_stone,EA,Ilat, &
      & spin,shape_spin,tableNN,shape_tableNN,masque,shape_masque,indexNN,shape_index,h_ext,my_lattice)
             E_new=sum(E_decompose)

             Spin(4:6,i_x_s,i_y_s,i_z_s,i_m_s)=S_old
             E_test=E_total+E_new-E_old

! position of the energy difference in energy space
   ! check if the energy is smaller that the minimum of energy
!             if ((E_test.lt.E_min).or.(E_test.gt.E_max)) cycle

   ! locate the position of E
             n_old=INT((E_total-E_min)/E_delta)+1
             n_new=INT((E_total+E_new-E_old-E_min)/E_delta)+1

   ! density of states
             lng_old=lng(n_old)
             lng_new=lng(n_new)
             delta_lng=lng_old-lng_new

   ! Do I accept the move or not???
            if (delta_lng.ge.1.0d2) then
                accept=.True.
            else
#ifdef CPP_MRG
                choice=mtprng_rand_real1(state)
#else
                CALL RANDOM_NUMBER(choice)
#endif
                if (choice.lt.1.0d-9) then
                  accept=.True.
                else
                  if (log(choice).lt.delta_lng) accept=.True.
                endif

            endif

            if (accept) then
                Spin(4:6,i_x_s,i_y_s,i_z_s,i_m_s)=S_new
                lng(n_new)=lng(n_new)+lnf
                count_E(n_new)=count_E(n_new)+1
                E_total=E_total+E_new-E_old
                accept=.False.
!                length=length-int(S_old(3))+int(S_new(3))
             else
                count_E(n_old)=count_E(n_old)+1
                lng(n_old)=lng(n_old)+lnf
             endif

!             vector((N_cell+length)/2+1)=vector((N_cell+length)/2+1)+1
! check the flatness of the histogramm
! The minimum of energy and the max must be taken out
             n_sweep=sum(count_E)
             mean_histo=real(n_sweep)/real(N_cell)
             min_histo=minval(count_E)

#ifdef CPP_DEBUG
             if (mod(i_step,100).eq.0) then
                 write(*,*) 'Check'
                 write(*,*) lng
                 write(*,*) count_E
                 write(*,*) min_histo,flatness*mean_histo,N_cell,n_sweep,mean_histo
             endif
#endif

             if ((mod(i_step,1000).eq.0).and.(real(min_histo).gt.flatness*mean_histo)) then
                 write(6,'(a,2x,I10,2x,a)') 'histogramm flat enough after',i_step,'steps'
                 write(6,'(a/)') 'reset the histogramm and the parameter f'
                 write(6,'(/a,I10/)') 'end of loop number', i_loop
                 exit
             endif

             if (i_step.eq.n_steps) then
                 write(6,'(a)') 'The simulation did not converge - Histogram not flat'
                 i_exit=.true.
                 exit
             endif
          enddo   ! end over the loop of the energy histogram

          write(6,'(/a,I10/)') "writing the results for loop", i_loop

          do i=1,n_histo
              write(7,'(3(E20.10E3,2x),I10)') (E_min+real(i-1)*E_delta)/real(N_cell),lng(i),real(count_E(i))/real(n_cell),count_E(i)
          enddo
          write(7,'(a)') ''

          if (i_exit) stop 'stoping the simulation'

! update of f and reset of the the energy histogram
          lnf=lnf/2.0d0
          count_E=0

      enddo

      close(7)

#ifdef CPP_DEBUG
      OPEN(666,FILE='test.dat',action='write',status='unknown',form='formatted')
      do i=1,N_cell+1
        write(666,'(I10)') vector(i)
      enddo
      close(666)
#endif

      write(6,'(/a)') "-------------------------"
      write(6,'(a/)') "writing the final results"

! renormalizing the density g(E)
!      lng=lng-log(real(N_cell))
      delta_lng=minval(lng)

      OPEN(7,FILE='lngE.dat',action='write',status='unknown',form='formatted')
      do i=1,n_histo
         write(7,'(2(E20.10E3,2x))') (E_min+real(i-1)*E_delta)/real(N_cell),lng(i)-delta_lng
      enddo
      close(7)

! calculate thermodynamical quantitites

      do i_T=1,n_Tsteps
          kT=(kTfin-kTini)/real(n_Tsteps-1)*real(i_T-1)+kTini
          temperatures(i_T)=kT/k_B
      enddo

      do i_T=1,n_Tsteps
          kT=temperatures(i_T)*k_B
          do i_histo=1,n_histo

                DE = (real(i_histo-1)*E_delta+E_min)
                DProb=lng(i_histo)-delta_lng-(DE+E_max)/kt

                Z(i_T)=Z(i_T)+exp(DProb)
                energy(i_T)=energy(i_T)+DE*exp(DProb)
                energy_sq(i_T)=energy_sq(i_T)+(DE)**2*exp(DProb)

          enddo
! renormalize the energy
          energy(i_T)=energy(i_T)/Z(i_T)
          energy_sq(i_T)=energy_sq(i_T)/Z(i_T)
          free_E(i_T)=-kT*log(Z(i_T))
          entropy(i_T)=k_B*(log(Z(i_T))+(energy(i_T)+E_max)/kT)
          heat_cap(i_T)=(energy_sq(i_T)-energy(i_T)**2)/kT**2*k_B
      enddo

      OPEN(7,FILE='EM.dat',action='write',status='unknown',form='formatted')
      write(7,'(a)') '#1:T   2:Z(T)   3:E(T)   4:F(B,T)   5:S(B,T)   6:C_V(B,T)'
      do i_T=1,n_Tsteps
          write(7,'(6(E20.10E3,2x))') temperatures(i_T),Z(i_T),energy(i_T),free_E(i_T),entropy(i_T),heat_cap(i_T)
      enddo
      close(7)

      write(6,'(a)') "done with entropic"

      end subroutine entropic
