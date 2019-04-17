      module m_parameters
      logical :: i_dispersion,i_sliptun,i_metropolis
      logical :: i_gneb,i_minimization
      logical :: i_qorien
      integer :: n_thousand
! Wang Landau sampling
      logical :: i_entropic
! size of the relaxation file
      integer :: n_sizerelax
! number of exchange layers as a parameter for intra and interlayer coupling
      integer :: param_N_Nneigh, param_N_Nneigh_il
!  Real*8, Parameter :: mu_B=0.0001156617d0
! exchange constants
      real(kind=8), dimension(:,:), allocatable :: J_ij
      real(kind=8), dimension(:), allocatable :: j_z
      real(kind=8), dimension(:), allocatable :: j_il
! Dzyaloshinsky-Moriya interaction
      real(kind=8), dimension(:,:), allocatable :: DM
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
! vectors for the DM TripleProduct
      real(kind=8), allocatable, Dimension(:,:,:) :: DM_vector
! frequence of plotting the graph for the spin dynamics
      integer :: gra_freq
! which programfiles should be written
      Logical :: Cor_log, Gra_log, gra_fft, gra_topo,dispersion
! Periodic boundary conditions and cool down or heat up algorithm
      Logical :: CalTheta,CalEnergy
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
subroutine inp_rw(io_simu,N_Nneigh,phase,Nei_z,Nei_il)
use m_parameters
use m_constants
use m_derived_types
use m_io_utils
use m_io_files_utils

#ifdef CPP_MPI
      use m_mpi_prop, only : irank,isize
      use mpi
#endif
      implicit none
!ccccccccccccccccccccccccccccccccccccccccccc
!In/out variable
real(kind=8), intent(inout) :: phase
integer, intent(inout) :: Nei_z,Nei_il
type(io_parameter), intent(inout) :: io_simu

! internal variables
integer  :: io_input
! local variables
      real(kind=8) :: DM_loc(12),J_loc(12),J_illoc(12)
!      real*8 :: n_spins_column1,n_spins_column2,n_spins_rows
!ccccccccccccccccccccccccccccccccccccccccccc
!in out variable
      integer, intent(inout) :: N_Nneigh 
! dummy variable
      integer :: fin,i,j,N_Nneigh_il
      character(len=100) :: str
      character(len=10) :: dummy
      logical :: exists,DM_userdef,asym
      character(len=4) :: tag
      real (kind=8) :: jz(12),jil(12)

      n_ghost=1
      nRepProc=1
      print_relax=.True.
      i_separate=.False.
      i_average=.False.
      i_ghost=.False.
      T_sweep=.False.

      i_optTset=.False.
      DM_userdef=.False.
      asym=.False.
      overrel=.False.
      underrel=.False.
      equi=.False.
      sphere=.True.
      dispersion=.False.
      i_sliptun=.False.
      i_qorien=.False.
      i_topohall=.False.
      CalTheta=.False.
      DM_loc=0.0d0
      J_loc=0.0d0
      J_illoc=0.0d0
      Jz=0.0d0
      Jil=0.0d0
      cone=2.0d0
      d_mu=1.0d0
      n_thousand=1000
      n_sizerelax=1000
      gra_freq=1
      N_temp=1
!      pi=acos(0.d0)*2.d0

io_input=open_file_read('input')

! mpi variables
call get_parameter(io_input,'input','ghost',i_ghost)
call get_parameter(io_input,'input','algo_mpi',i_ghost)


call get_parameter(io_input,'input','nRepProc',nRepProc)
call get_parameter(io_input,'input','n_sizerelax',n_sizerelax)
call get_parameter(io_input,'input','n_relaxation',n_thousand)

! io variables
call get_parameter(io_input,'input','gra_fft',io_simu%io_fft_Xstruct)
call get_parameter(io_input,'input','gra_topo',io_simu%io_topo)
call get_parameter(io_input,'input','qorien',io_simu%io_qorien)
call get_parameter(io_input,'input','warnings',io_simu%io_warning)
call get_parameter(io_input,'input','dispersion',io_simu%io_dispersion)
call get_parameter(io_input,'input','gra_log',io_simu%io_Xstruct)
call get_parameter(io_input,'input','gra_freq',io_simu%io_frequency)
call get_parameter(io_input,'input','SPSTM-only',io_simu%io_spstmonly)
call get_parameter(io_input,'input','SPSTM-image',io_simu%io_spstmL)

!do i=1,12
!  call convert()
!  call get_parameter(io_input,'input','SPSTM-image',io_simu%io_spstmL)
!enddo

      rewind(io_input)
      do
      read (io_input,'(a)',iostat=fin) str
        if (fin /= 0) exit
        str= trim(adjustl(str))
        if (len_trim(str)==0) cycle
        if (str(1:1) == '#' ) cycle 
!cccc We start to read the input


         if ( str(1:8) == 'algo_mpi') then
           backspace(io_input)
           read(io_input,*) dummy, tag
           if (tag == "sepa") then
            i_separate=.True.
           elseif (tag == "aver") then
            i_average=.True.
           elseif (tag == "para") then
            i_separate=.True.
            backspace(io_input)
            read(io_input,*) dummy, dummy, N_temp, i_optTset, T_relax_temp
           elseif (tag == "none") then
            i_separate=.False.
            i_average=.False.
           else

           endif
         endif

         if ( str(1:8) == 'sampling') then
           backspace(io_input)
           read(io_input,*) dummy, tag
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
! read the constant from the inp file

        if ( str(1:8) == 'CalTheta') then
           backspace(io_input)
           read(io_input,*) dummy, CalTheta
          endif
        if ( str(1:9) == 'CalEnergy') then
           backspace(io_input)
           read(io_input,*) dummy, CalEnergy
          endif

        if ( str(1:4) == 'cone') then
           backspace(io_input)
           read(io_input,*) dummy, cone
          endif
        if ( str(1:4) == 'd_mu') then
           backspace(io_input)
           read(io_input,*) dummy, d_mu
          endif

!           kTini=kTini*k_B

!!! Interlayer coupling
        if ( str(1:6) == 'J_1il') then
           i_sliptun=.True.
           backspace(io_input)
           read(io_input,*) dummy, J_illoc(1)
           if (dabs(J_illoc(1)).gt.eps) then
           N_Nneigh_il=1
           endif
          endif
        if ( str(1:6) == 'J_2il') then
           i_sliptun=.True.
           backspace(io_input)
           read(io_input,*) dummy, J_illoc(2)
           if (dabs(J_illoc(2)).gt.eps) then
           N_Nneigh_il=2
           endif
          endif
        if ( str(1:6) == 'J_3il') then
           i_sliptun=.True.
           backspace(io_input)
           read(io_input,*) dummy, J_illoc(3)
           if (dabs(J_illoc(3)).gt.eps) then
           N_Nneigh_il=3
           endif
          endif
        if ( str(1:6) == 'J_4il') then
           i_sliptun=.True.
           backspace(io_input)
           read(io_input,*) dummy, J_illoc(4)
           if (dabs(J_illoc(4)).gt.eps) then
           N_Nneigh_il=4
           endif
          endif
        if ( str(1:6) == 'J_5il') then
           i_sliptun=.True.
           backspace(io_input)
           read(io_input,*) dummy, J_illoc(5)
           if (dabs(J_illoc(5)).gt.eps) then
           N_Nneigh_il=5
           endif
          endif
        if ( str(1:6) == 'J_6il') then
           i_sliptun=.True.
           backspace(io_input)
           read(io_input,*) dummy, J_illoc(6)
           if (dabs(J_illoc(6)).gt.eps) then
           N_Nneigh_il=6
           endif
          endif
        if ( str(1:6) == 'J_7il') then
           i_sliptun=.True.
           backspace(io_input)
           read(io_input,*) dummy, J_illoc(7)
           if (dabs(J_illoc(7)).gt.eps) then
           N_Nneigh_il=7
           endif
          endif
        if ( str(1:6) == 'J_8il') then
           i_sliptun=.True.
           backspace(io_input)
           read(io_input,*) dummy, J_illoc(8)
           if (dabs(J_illoc(8)).gt.eps) then
           N_Nneigh_il=8
           endif
          endif
        if ( str(1:6) == 'J_9il') then
           i_sliptun=.True.
           backspace(io_input)
           read(io_input,*) dummy, J_illoc(9)
           if (dabs(J_illoc(9)).gt.eps) then
           N_Nneigh_il=9
           endif
          endif
        if ( str(1:6) == 'J_10il') then
           i_sliptun=.True.
           backspace(io_input)
           read(io_input,*) dummy, J_illoc(10)
           if (dabs(J_illoc(10)).gt.eps) then
           N_Nneigh_il=10
           endif
          endif

!!! save the highest order of intralayer J, that is unequal zero in "inp" as a parameter
       param_N_Nneigh=N_Nneigh


        if ( str(1:14) == 'Total_MC_Steps') then
           backspace(io_input)
           read(io_input,*) dummy, Total_MC_Steps
          endif
        if ( str(1:8) == 'n_Tsteps') then
           backspace(io_input)
           read(io_input,*) dummy, n_Tsteps
          endif
        if ( str(1:7) == 'T_relax') then
           backspace(io_input)
           read(io_input,*) dummy, T_relax
          endif
        if ( str(1:6) == 'T_auto') then
           backspace(io_input)
           read(io_input,*) dummy, T_auto
          endif
        if ( str(1:7) == 'Cor_log') then
           backspace(io_input)
           read(io_input,*) dummy, Cor_log
          endif

        if ( str(1:10) == 'relaxation') then
           backspace(io_input)
           read(io_input,*) dummy, tag
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

call close_file('input',io_input)

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
      endif

    ! change by SvM: allocate DMI parameters
    if (count(dabs(DM_loc).gt.1.0d-8).ne.0) then
        j=0
        do i=1,12
            if (dabs(DM_loc(i)).gt.1.0d-8) j=i
        enddo
        allocate(DM(j,1))
        DM(1:j,1)=DM_loc(1:j)
    else
        allocate(DM(1,1))
        DM=0.0d0
    endif

      inquire (file='zdir.in',exist=exists)
      if (exists) then
       call rw_zdir(phase,jz,Nei_z,jil,Nei_il)
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


      if (n_sizerelax.gt.n_thousand) then
       write(6,'(a)') 'please choose n_sizerelax < n_relaxation'
      endif

#ifdef CPP_DEBUG
      write(6,*) "nb neighbor", N_Nneigh
      write(6,*) "J_ij", J_ij
      write(6,*) "J_z", Jz
#endif

! put local variables into the types for transfer

      end subroutine inp_rw

