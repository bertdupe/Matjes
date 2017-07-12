      module m_parameters
      logical :: i_dispersion,i_sliptun,i_metropolis
      logical :: i_biq,i_four,i_DM,i_stone,i_gneb,i_minimization
      logical :: i_qorien
      integer :: n_thousand
! Wang Landau sampling
      logical :: i_entropic
! size of the relaxation file
      integer :: n_sizerelax
! number of exchange layers as a parameter for intra and interlayer coupling
      integer :: param_N_Nneigh, param_N_Nneigh_il
! check for dipole dipole contribution
      logical :: i_dip
! Temperature of the heat bath in the units [kT/J]
     real(kind=8) :: kT, kTfin, kTini
!  Real*8, Parameter :: mu_B=0.0001156617d0
!  external magnetic field
      real(kind=8), Dimension(1:3) :: H_ext
! exchange constants
      real(kind=8), dimension(:,:), allocatable :: J_ij
      real(kind=8), dimension(:), allocatable :: j_z
      real(kind=8), dimension(:), allocatable :: j_il
! Dzyaloshinsky-Moriya interaction
      real(kind=8), dimension(:,:), allocatable :: DM
! easy axis
      real(kind=8) :: EA(3)
! Go in the spmstm program of Tobias
      Logical :: spstmL, spstmonly
! sampling system
      logical :: sphere,equi,dynamic
! relaxation of the system
      logical :: overrel, underrel
! cone angle
      real(kind=8) :: cone,d_mu
! Total number of the Monte Carlo Steps
      INTEGER :: Total_MC_Steps
! Total number of temperature steps
      Integer :: n_Tsteps
! Relaxation times in case of the last step or all others
      Integer :: T_relax
! Autocorrelation Time
      Integer :: T_auto
! number of figures where the spinconfiguration is shown
      Integer :: n_spingra
!stoner parameter
      real(kind=8) :: Ist
! biquadratic spin coupling constant
      real(kind=8) :: J_B
! vectors for the DM TripleProduct
      real(kind=8), allocatable, Dimension(:,:,:) :: DM_vector
! constants in the sum
      real(kind=8) :: c_Ji,c_DM,c_JB,c_Ki,c_ani
! Anisotropie Vector
      real(kind=8), Dimension(1:3) :: D_ani
! four spin interaction
      real(kind=8) :: K_1
! frequence of plotting the graph for the spin dynamics
      integer :: gra_freq
! which programfiles should be written
      Logical :: Cor_log, Gra_log, gra_fft, gra_topo,dispersion
! Periodic boundary conditions and cool down or heat up algorithm
      Logical :: Periodic_log(3),i_cool,CalTheta,CalEnergy
! calculate the topological Hall effect
      logical ::  i_topohall
! epsilon
      real(kind=8), parameter :: eps=1.0d-8
! ising model
      logical :: ising
! print a lots of warnings if necessary
      logical :: i_print_W
      logical :: i_separate,i_average,i_ghost,i_optTset,print_relax
      logical :: T_sweep
      integer :: N_temp,T_relax_temp
      integer :: n_ghost,nRepProc
#ifdef CPP_MPI
      integer :: cell_irank
#endif
      end module

!!!!!!!!!!!!!!!!!!!!!!!!!!
! This routine reads the input file and setup and every variables given in the inp file
!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine inp_rw(my_simu,io_simu,N_Nneigh,phase,Nei_z,Nei_il)
      use m_parameters
      use m_rw_lattice, only : net,dim_lat
      use m_constants
      use m_derived_types
#ifdef CPP_MPI
      use m_mpi_prop, only : irank,isize
      use mpi
#endif
      implicit none
!ccccccccccccccccccccccccccccccccccccccccccc
!In/out variable
      real(kind=8), intent(inout) :: phase
      integer, intent(inout) :: Nei_z,Nei_il
      type(io_parameter), intent(out) :: io_simu
      type(type_simu), intent(out) :: my_simu

! local variables
      real(kind=8) :: angle,v(3),vort,qeuler
      real(kind=8) :: alph30,betha,SQRT1d2
      real(kind=8) :: DM_loc(3),J_loc(12),J_illoc(12)
!      real*8 :: n_spins_column1,n_spins_column2,n_spins_rows
!ccccccccccccccccccccccccccccccccccccccccccc
!in out variable
      integer, intent(inout) :: N_Nneigh 
! dummy variable
      integer, parameter  :: io=1
      integer :: fin,i,j,ierr,N_Nneigh_il
      character(len=100) :: str
      character(len=10) :: dummy
      logical :: exists,DM_userdef,asym,topoonly,i_exi
      character(len=1) :: toto
      character(len=4) :: tag
      real (kind=8) :: jz(12),jil(12)

! type of the simulation
      my_simu%i_paratemp=.False.
      my_simu%i_dynamic=.False.
      my_simu%i_metropolis=.False.
      my_simu%i_gneb=.False.
      my_simu%i_minimization=.False.
      my_simu%i_entropic=.False.
      my_simu%i_r_texture=.False.

! io_of the simulation

      n_ghost=1
      nRepProc=1
      print_relax=.True.
      i_separate=.False.
      i_average=.False.
      i_ghost=.False.
      T_sweep=.False.

      i_optTset=.False.
      i_DM=.False.
      ising=.False.
      DM_userdef=.False.
      asym=.False.
      overrel=.False.
      underrel=.False.
      equi=.False.
      sphere=.True.
      dispersion=.False.
      i_biq=.False.
      I_four=.False.
      i_stone=.False.
      i_sliptun=.False.
      i_gneb=.False.
      i_qorien=.False.
      i_topohall=.False.
      i_print_W=.False.
      CalTheta=.False.
      Gra_log=.False.
      DM_loc=0.0d0
      J_loc=0.0d0
      J_illoc=0.0d0
      Jz=0.0d0
      Jil=0.0d0
      cone=2.0d0
      d_mu=1.0d0
      J_B=0.0d0
      kT=0.0d0
      kTfin=10.0d0
      kTini=1.0d0
      kT=kTini
      n_thousand=1000
      n_sizerelax=1000
      EA=(/0,0,1/)
      gra_freq=1
      N_temp=1


      c_Ji=-1.0d0
      c_DM=-1.0d0
      c_JB=-1.0d0
      c_Ki=-1.0d0
      c_ani=1.0d0
!      pi=acos(0.d0)*2.d0

! open the input
      inquire (file='inp',exist=exists)
      if (.not. exists) then
      write(6,*) 'no input file'
      STOP
      endif

      open (io,file='inp',form='formatted',status='old',action='read')

      rewind(io)
      do
      read (io,'(a)',iostat=fin) str
        if (fin /= 0) exit
        str= trim(adjustl(str))
        if (len_trim(str)==0) cycle

        if (str(1:1) == '#' ) cycle 
!cccc We start to read the input

         if ( str(1:10) == 'simulation') then
           backspace(io)
           read(io,*) dummy, dummy
           select case (adjustl(dummy))
            case ("Montecarlo")
              my_simu%i_metropolis=.True.
            case ("spindynami")
              my_simu%i_dynamic=.True.
            case ("entropic")
              my_simu%i_entropic=.True.
            case ("GNEB")
              my_simu%i_gneb=.True.
            case ("parallelte")
              my_simu%i_paratemp=.True.
            case ("minimizati")
              my_simu%i_minimization=.True.
            case default
               STOP 'select a simulation type'
           end select
         endif

         if ( str(1:5) == 'ghost') then
           backspace(io)
           read(io,*) dummy, i_ghost, n_ghost
         endif
         if ( str(1:8) == 'nRepProc') then
           backspace(io)
           read(io,*) dummy, nRepProc
         endif
         if ( str(1:8) == 'algo_mpi') then
           backspace(io)
           read(io,*) dummy, tag
           if (tag == "sepa") then
            i_separate=.True.
           elseif (tag == "aver") then
            i_average=.True.
           elseif (tag == "para") then
            i_separate=.True.
            backspace(io)
            read(io,*) dummy, dummy, N_temp, i_optTset, T_relax_temp
           elseif (tag == "none") then
            i_separate=.False.
            i_average=.False.
           else
#ifdef CPP_MPI
            if (irank.eq.0) then
#endif
            write(6,'(a)') 'cannot read algo_mpi algo'
#ifdef CPP_MPI
            call mpi_finalize(ierr)
            endif
#else
            stop
#endif
           endif
         endif
        if (my_simu%i_paratemp.or.(str(1:11) == 'print_relax')) then
         backspace(io)
         read(io,*) dummy, print_relax
        endif
        if ( str(1:11) == 'n_sizerelax') then
           backspace(io)
           read(io,*) dummy, n_sizerelax
          endif
        if ( str(1:5) == 'ising') ising=.True.
         if ( str(1:7) == 'T_sweep') then
           backspace(io)
           read(io,*) dummy,T_sweep
         endif
         if ( str(1:8) == 'sampling') then
           backspace(io)
           read(io,*) dummy, tag 
           if (tag == "equi") then
            equi=.True.
            sphere=.False.
           elseif (tag == "sphe") then
            equi=.False.
            sphere=.True.
           else
            write(*,*) "choose a sampling"
            STOP
           endif
         endif
        if ( str(1:6) == 'stoner') then
           backspace(io)
           i_stone=.True.
           read(io,*) dummy, Ist
          endif
! read the constant from the inp file
        if ( str(1:12) == 'n_relaxation') then
           backspace(io)
           read(io,*) dummy, n_thousand
          endif
        if ( str(1:4) == 'c_Ji') then
           backspace(io)
           read(io,*) dummy, c_Ji
          endif
        if ( str(1:4) == 'c_DM') then
           backspace(io)
           read(io,*) dummy, c_DM
          endif
        if ( str(1:4) == 'c_JB') then
           backspace(io)
           read(io,*) dummy, c_JB
          endif
        if ( str(1:4) == 'c_Ki') then
           backspace(io)
           read(io,*) dummy, c_Ki
          endif
        if ( str(1:5) == 'c_ani') then
           backspace(io)
           read(io,*) dummy, c_ani
          endif
        if ( str(1:7) == 'gra_fft') then
           backspace(io)
           read(io,*) dummy, gra_fft
          endif
        if ( str(1:8) == 'CalTheta') then
           backspace(io)
           read(io,*) dummy, CalTheta
          endif
        if ( str(1:9) == 'CalEnergy') then
           backspace(io)
           read(io,*) dummy, CalEnergy
          endif
        if ( str(1:8) == 'gra_topo') then
           backspace(io)
           read(io,*) dummy, gra_topo
          endif
        if ( str(1:6) == 'qorien') then
           backspace(io)
           read(io,*) dummy, i_qorien
          endif
        if ( str(1:4) == 'cone') then
           backspace(io)
           read(io,*) dummy, cone
          endif
        if ( str(1:4) == 'd_mu') then
           backspace(io)
           read(io,*) dummy, d_mu
          endif
        if ( str(1:5) == 'H_ext') then
           backspace(io)
           read(io,*) dummy, (H_ext(i),i=1,3)
          endif
        if ( str(1:2) == 'EA') then
           backspace(io)
           read(io,*) dummy, (EA(i),i=1,3)
          endif
        if ( str(1:4) == 'Tini') then
           backspace(io)
           read(io,*) dummy, kTini
           kTini=kTini*k_B
           kT=kTini
          endif
        if ( str(1:4) == 'Tfin') then
           backspace(io)
           read(io,*) dummy, kTfin
           kTfin=kTfin*k_B
          endif

!!! Interlayer coupling
        if ( str(1:6) == 'J_1il') then
           i_sliptun=.True.
           backspace(io)
           read(io,*) dummy, J_illoc(1)
           if (dabs(J_illoc(1)).gt.eps) then
           N_Nneigh_il=1
           endif
          endif
        if ( str(1:6) == 'J_2il') then
           i_sliptun=.True.
           backspace(io)
           read(io,*) dummy, J_illoc(2)
           if (dabs(J_illoc(2)).gt.eps) then
           N_Nneigh_il=2
           endif
          endif
        if ( str(1:6) == 'J_3il') then
           i_sliptun=.True.
           backspace(io)
           read(io,*) dummy, J_illoc(3)
           if (dabs(J_illoc(3)).gt.eps) then
           N_Nneigh_il=3
           endif
          endif
        if ( str(1:6) == 'J_4il') then
           i_sliptun=.True.
           backspace(io)
           read(io,*) dummy, J_illoc(4)
           if (dabs(J_illoc(4)).gt.eps) then
           N_Nneigh_il=4
           endif
          endif
        if ( str(1:6) == 'J_5il') then
           i_sliptun=.True.
           backspace(io)
           read(io,*) dummy, J_illoc(5)
           if (dabs(J_illoc(5)).gt.eps) then
           N_Nneigh_il=5
           endif
          endif
        if ( str(1:6) == 'J_6il') then
           i_sliptun=.True.
           backspace(io)
           read(io,*) dummy, J_illoc(6)
           if (dabs(J_illoc(6)).gt.eps) then
           N_Nneigh_il=6
           endif
          endif
        if ( str(1:6) == 'J_7il') then
           i_sliptun=.True.
           backspace(io)
           read(io,*) dummy, J_illoc(7)
           if (dabs(J_illoc(7)).gt.eps) then
           N_Nneigh_il=7
           endif
          endif
        if ( str(1:6) == 'J_8il') then
           i_sliptun=.True.
           backspace(io)
           read(io,*) dummy, J_illoc(8)
           if (dabs(J_illoc(8)).gt.eps) then
           N_Nneigh_il=8
           endif
          endif
        if ( str(1:6) == 'J_9il') then
           i_sliptun=.True.
           backspace(io)
           read(io,*) dummy, J_illoc(9)
           if (dabs(J_illoc(9)).gt.eps) then
           N_Nneigh_il=9
           endif
          endif
        if ( str(1:6) == 'J_10il') then
           i_sliptun=.True.
           backspace(io)
           read(io,*) dummy, J_illoc(10)
           if (dabs(J_illoc(10)).gt.eps) then
           N_Nneigh_il=10
           endif
          endif

!!! save the highest order of interlayer J, that is unequal zero in "inp" as a parameter
       param_N_Nneigh_il=N_Nneigh_il

!! Intralayer coupling
        if ( str(1:4) == 'J_1') then
           backspace(io)
           read(io,*) dummy, J_loc(1)
           if (dabs(J_loc(1)).gt.eps) then
           N_Nneigh=1
           endif
          endif
        if ( str(1:4) == 'J_2') then
           backspace(io)
           read(io,*) dummy, J_loc(2)
           if (dabs(J_loc(2)).gt.eps) then
           N_Nneigh=2
           endif
          endif
        if ( str(1:4) == 'J_3') then
           backspace(io)
           read(io,*) dummy, J_loc(3)
           if (dabs(J_loc(3)).gt.eps) then
           N_Nneigh=3
           endif
          endif
        if ( str(1:4) == 'J_4') then
           backspace(io)
           read(io,*) dummy, J_loc(4)
           if (dabs(J_loc(4)).gt.eps) then
           N_Nneigh=4
           endif
          endif
        if ( str(1:4) == 'J_5') then
           backspace(io)
           read(io,*) dummy, J_loc(5)
           if (dabs(J_loc(5)).gt.eps) then
           N_Nneigh=5
           endif
          endif
        if ( str(1:4) == 'J_6') then
           backspace(io)
           read(io,*) dummy, J_loc(6)
           if (dabs(J_loc(6)).gt.eps) then
           N_Nneigh=6
           endif
          endif
        if ( str(1:4) == 'J_7') then
           backspace(io)
           read(io,*) dummy, J_loc(7)
           if (dabs(J_loc(7)).gt.eps) then
           N_Nneigh=7
           endif
           endif
        if ( str(1:4) == 'J_8') then
           backspace(io)
           read(io,*) dummy, J_loc(8)
           if (dabs(J_loc(8)).gt.eps) then
           N_Nneigh=8
           endif
           endif
        if ( str(1:4) == 'J_9') then
           backspace(io)
           read(io,*) dummy, J_loc(9)
           if (dabs(J_loc(9)).gt.eps) then
           N_Nneigh=9
           endif
           endif
        if ( str(1:4) == 'J_10') then
           backspace(io)
           read(io,*) dummy, J_loc(10)
           if (dabs(J_loc(10)).gt.eps) then
           N_Nneigh=10
           endif
           endif
        if ( str(1:4) == 'J_11') then
           backspace(io)
           read(io,*) dummy, J_loc(11)
           if (dabs(J_loc(11)).gt.eps) then
           N_Nneigh=11
           endif
           endif
        if ( str(1:4) == 'J_12') then
           backspace(io)
           read(io,*) dummy, J_loc(12)
           if (dabs(J_loc(12)).gt.eps) then
           N_Nneigh=12
           endif

           endif

!!! save the highest order of intralayer J, that is unequal zero in "inp" as a parameter
       param_N_Nneigh=N_Nneigh

       if ( str(1:3) == 'J_B') then
           backspace(io)
           read(io,*) dummy, J_B
           if (dabs(J_B).gt.1.0d-8) i_biq=.True.
          endif
        if ( str(1:6) == 'dipdip') then
           backspace(io)
           read(io,*) dummy, i_dip
          endif
        if ( str(1:5) == 'D_ani') then
           backspace(io)
           read(io,*) dummy, (D_ani(i),i=1,3)
          endif
        if ( str(1:3) == 'K_1') then
           backspace(io)
           read(io,*) dummy, K_1
           if (dabs(K_1).gt.1.0d-8) i_four=.True.
          endif
        if ( str(1:8) == 'warnings') then
           backspace(io)
           read(io,*) dummy, i_print_W
          endif

        if (( str(1:2) == 'DM').and.(.not.i_dm)) then
           backspace(io)
          if ( str(1:4) == 'DM_1') then
           i_DM=.True.
           read(io,*) dummy, DM_loc(1), DM_loc(2), DM_loc(3)
           j=maxloc(abs(DM_loc),1)
           read(io,*) dummy, DM_loc(1), DM_loc(2), DM_loc(3)
           i=maxloc(abs(DM_loc),1)
           j=maxval((/i,j/))
           allocate(DM(j,2))
           DM(1:j,2)=DM_loc(1:j)
           backspace(io)
           backspace(io)
           read(io,*) dummy, DM_loc(1), DM_loc(2), DM_loc(3)
           DM(1:j,1)=DM_loc(1:j)
          else
           read(io,*) dummy, DM_loc(1), DM_loc(2), DM_loc(3)
           if (count(dabs(DM_loc).gt.1.0d-8).ne.0) then
            i_DM=.True.
            j=0
            do i=1,3
             if (dabs(DM_loc(i)).gt.1.0d-8) j=i
            enddo
            allocate(DM(j,1))
            DM(1:j,1)=DM_loc(1:j)
            else
            allocate(DM(1,1))
            DM=0.0d0
           endif
          endif
          endif

        if ( str(1:9) == 'DM-vector') then
           backspace(io)
           read(io,*) dummy, v(1), v(2), v(3)
           DM_userdef=.True.
           write(*,*) 'DM vector defined by user'
           write(*,*) 'I hope you know what you are doing'
          endif
        if ( str(1:14) == 'Total_MC_Steps') then
           backspace(io)
           read(io,*) dummy, Total_MC_Steps
          endif
        if ( str(1:8) == 'n_Tsteps') then
           backspace(io)
           read(io,*) dummy, n_Tsteps
          endif
        if ( str(1:7) == 'T_relax') then
           backspace(io)
           read(io,*) dummy, T_relax
          endif
        if ( str(1:6) == 'T_auto') then
           backspace(io)
           read(io,*) dummy, T_auto
          endif
        if ( str(1:9) == 'n_spingra') then
           backspace(io)
           read(io,*) dummy, n_spingra
          endif
        if ( str(1:7) == 'Cor_log') then
           backspace(io)
           read(io,*) dummy, Cor_log
          endif
        if ( str(1:7) == 'Gra_log') then
           backspace(io)
           read(io,*) dummy, Gra_log, gra_freq
          endif
        if ( str(1:12) == 'Periodic_log') then
           backspace(io)
           read(io,*) dummy, Periodic_log(1), Periodic_log(2), Periodic_log(3)
          endif
        if ( str(1:10) == 'dispersion') then
           backspace(io)
           read(io,*) dummy, dispersion
          endif
        if ( str(1:10) == 'relaxation') then
           backspace(io)
           read(io,*) dummy, tag
           if (tag == "over") then
           overrel=.True.
           elseif (tag == "unde") then
           underrel=.True.
           else
           underrel=.False.
           overrel=.False.
           endif
          endif
      end do
      close(io)

      if (count(dabs(J_loc).gt.1.0d-8).ne.0) then
       j=0
       do i=1,12
        if (dabs(J_loc(i)).gt.1.0d-8) j=i
       enddo

       if (i_sliptun) then
        allocate(J_ij(j,2))
        J_ij(:,1)=J_loc(1:j)
        do i=1,j
         if (abs(J_illoc(i)).gt.0.0d0) then
          J_ij(i,2)=J_illoc(i)
          else
          J_ij(i,2)=J_loc(i)
         endif
        enddo
        else
        allocate(J_ij(j,1))
        J_ij(:,1)=J_loc(1:j)
       endif
      else
       write(6,*) "no J_ij have been found"
       stop
      endif
      i_cool=.false.
      if (kTini.gt.kTfin) i_cool=.true.

      inquire (file='zdir.in',exist=exists)
      if (exists) then
       call rw_zdir(phase,jz,Nei_z,jil,Nei_il,dim_lat)
       if (Nei_il.ne.0) then
        allocate(j_il(Nei_il))
        j_il=jil(1:Nei_il)
       else
        allocate(j_il(1))
        j_il=0.0d0
       endif
       if (Nei_z.ne.0) then
        allocate(j_z(Nei_z))
        j_z=jz(1:Nei_z)
       else
        allocate(j_z(1))
        j_z=0.0d0
       endif
       else
       allocate(j_z(1))
       allocate(j_il(1))
       j_z=0.0d0
       j_il=0.0d0
      endif


! check for electric field
      call rw_efield(dim_lat,net)
!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (n_sizerelax.gt.n_thousand) then
       write(6,'(a)') 'please choose n_sizerelax < n_relaxation'
      endif
#ifdef CPP_MPI
      if ((irank.eq.0).AND.(isize.lt.9)) then
       call SignatureFile(1,len_trim('param-'//char(48+irank)//'.dat'), &
        'param-'//char(48+irank)//'.dat',len_trim('asis'),'asis')
        elseif (isize.lt.9) then
       call SignatureFile(1,len_trim('param-'//char(48+irank)//'.dat'), &
        'param-'//char(48+irank)//'.dat',len_trim('asis'),'asis')
      endif
#else
       call SignatureFile(1,len_trim('param.dat'), &
        'param.dat',len_trim('asis'),'asis')
#endif

#ifdef CPP_DEBUG
      write(6,*) "J_ij", J_ij
      write(6,*) "J_z", Jz
#endif

! put local variables into the types for transfer

      end subroutine inp_rw

