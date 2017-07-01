!
! ===============================================================
      SUBROUTINE Relaxation(N_cell,n_system,kT,state,E_total,E,magnetization,qeulerp,qeulerm,vortex, &
    &  n_relaxation,n_sizerelax,T_relax,acc,rate,tries,cone,print_relax,h_ext,EA, &
    &  spin,shape_spin,tableNN,shape_tableNN,masque,shape_masque,indexNN,shape_index, &
    &  i_biq,i_dip,i_DM,i_four,i_stone,ising,i_print_W,equi,overrel,sphere,underrel,n_world,i_ghost)
      use mtprng
      use m_Corre
      use m_constants, only : k_b
      use m_topocharge_all
      use m_store_relaxation
#ifdef CPP_MPI
      use m_mpi_prop, only : isize,irank_working,MPI_COMM,MPI_COMM_BOX,irank_box
#endif
      Implicit none
! part of the interface
      interface
        subroutine MCsteps(state,E_total,E,Magnetization,kt,acc,rate,tries,cone,n_system, &
                & spin,shape_spin,tableNN,shape_tableNN,masque,shape_masque,indexNN,shape_index,h_ext,EA, &
                & i_biq,i_dip,i_DM,i_four,i_stone,ising,i_print_W,equi,overrel,sphere,underrel,n_world)
           use mtprng
           integer, intent(in) :: shape_index(2),shape_spin(5),shape_tableNN(6),shape_masque(4),n_system
           integer, intent(in) :: tableNN(shape_tableNN(1),shape_tableNN(2),shape_tableNN(3),shape_tableNN(4),shape_tableNN(5),shape_tableNN(6))
           integer, intent(in) :: masque(shape_masque(1),shape_masque(2),shape_masque(3),shape_masque(4))
           integer, intent(in) :: indexNN(shape_index(1),shape_index(2)),n_world
           logical, intent(in) :: i_biq,i_dip,i_DM,i_four,i_stone,ising,i_print_W,equi,overrel,sphere,underrel
           real(kind=8), intent(in) :: kt,h_ext(3),EA(3)
           real(kind=8), intent(inout) :: spin(shape_spin(1),shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5))
           type(mtprng_state), intent(inout) :: state
           real(kind=8), intent(inout) :: E_total,Magnetization(3),E(8),acc,rate,cone,tries
        end subroutine
      end interface
! input
      integer, intent(in) :: shape_index(2),shape_spin(5),shape_tableNN(6),shape_masque(4)
      real(kind=8), intent(inout) :: spin(shape_spin(1),shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5))
      type(mtprng_state), intent(inout) :: state
      real(kind=8), intent(inout) :: qeulerp,qeulerm,vortex(3),cone,acc,rate,tries
      real(kind=8), intent(inout) :: E_total,magnetization(3),E(8)
      real(kind=8), intent(in) :: kT,h_ext(3),EA(3)
      integer, intent(in) :: tableNN(shape_tableNN(1),shape_tableNN(2),shape_tableNN(3),shape_tableNN(4),shape_tableNN(5),shape_tableNN(6))
      integer, intent(in) :: masque(shape_masque(1),shape_masque(2),shape_masque(3),shape_masque(4))
      integer, intent(in) :: indexNN(shape_index(1),shape_index(2)),n_world
      integer, intent(in) :: n_relaxation,T_relax,N_cell,n_sizerelax,n_system
      logical, intent(in) :: print_relax,i_biq,i_dip,i_DM,i_four,i_stone,ising,i_print_W,equi,overrel,sphere,underrel,i_ghost
!    local variables. Used only in MPI setup
      real(kind=8) :: E_int(8),M_int(3),qp_int,qm_int
! a big table
      real(kind=8) :: Relax(18,n_sizerelax)
!     Slope Index
      integer :: i_relaxation,i_MC
! dummy
      integer :: i,j,i_pos,n_w_step,ierr
      character(len=30) :: fname,toto

#ifdef CPP_MPI

      include 'mpif.h'

      if (irank_working.eq.0) write(6,'(/,a,/)') "starting relaxation"
#else
      write(6,'(/,a,/)') "starting relaxation"
#endif

      Relax=0.0d0
      n_w_step=n_relaxation/n_sizerelax
!     Monte Carlo steps for thermal equilibrium
!     The first time it will take longer
!     -----------------------------------------------------------------

      do i_relaxation=1,n_relaxation
!         T_relax_1 is probably larger then T_relax, this is because
!         the last step might take more time to relax than an
!         the step to an unordered structure

            !Relaxation of the System
            Do i_MC=1,T_relax*N_cell
                Call MCStep(state,E_total,E,Magnetization,kt,acc,rate,tries,cone,n_system, &
                & spin,shape_spin,tableNN,shape_tableNN,masque,shape_masque,indexNN,shape_index,h_ext,EA, &
                & i_biq,i_dip,i_DM,i_four,i_stone,ising,i_print_W,equi,overrel,sphere,underrel,n_world)
            enddo

            !In case T_relax set to zero at least one MCstep is done
            Call MCStep(state,E_total,E,Magnetization,kt,acc,rate,tries,cone,n_system, &
            & spin,shape_spin,tableNN,shape_tableNN,masque,shape_masque,indexNN,shape_index,h_ext,EA, &
            & i_biq,i_dip,i_DM,i_four,i_stone,ising,i_print_W,equi,overrel,sphere,underrel,n_world)

! calculate the topocharge
            call topo(spin,shape_spin,masque,shape_masque,qeulerp,qeulerm)

            ! Write the Equilibrium files

            if (print_relax) then
#ifdef CPP_MPI
                if (i_ghost) then

                    E_int=E
                    M_int=Magnetization
                    qp_int=qeulerp
                    qm_int=qeulerm
                    qeulerm=0.0d0
                    qeulerp=0.0d0
                    Magnetization=0.0d0
                    E=0.0d0

                    call mpi_reduce(E,E_int,8,MPI_REAL8,MPI_SUM,0,MPI_COMM_BOX,ierr)
                    call mpi_reduce(Magnetization,M_int,3,MPI_REAL8,MPI_SUM,0,MPI_COMM_BOX,ierr)
                    call mpi_reduce(qeulerp,qp_int,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_BOX,ierr)
                    call mpi_reduce(qeulerm,qm_int,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_BOX,ierr)

                endif
#endif

                call store_relaxation(Relax,n_relaxation,i_relaxation,dble(i_relaxation), &
             &  sum(E)/dble(N_cell),E,dble(N_cell),kt,Magnetization,rate,cone,qeulerp,qeulerm)

            endif

      enddo   ! enddo over the relaxation loop
!!!!!!!!!!!!!***************************!!!!!!
!!!!!!!!!!!!!***************************!!!!!!

    write(6,'(/,a,f8.4,2x,a,/)') 'System is relaxed for T= ',kT/k_B,'Kelvin'

!print the Equilibrium files
      if (print_relax) then

! numbering of the files
          write(fname,'(f8.4)') kT/k_B
          toto=trim(adjustl(fname))

          write(fname,'(a,14a,a)')'Equilibriumphi-',(toto(i:i),i=1,len_trim(toto)),'.dat'
          OPEN(8,FILE=fname,status='unknown')

          Write(8,'(a)') "#   1:n_MC  2:Etot  3:T  4:M  5:rate  6:cone  7:Q  8:Exch  9:Zeeman  10:Ani &
   &         11:4S  12:DM  13:biq  14:dip  15:stoner  16:Chi(E)(t,T)  17:Chi(M)(t,T)  18:Chi(T)(t)"

          Relax(16,:)=correlation(Relax(2,:),n_sizerelax)
          Relax(17,:)=correlation(Relax(4,:),n_sizerelax)
          Relax(18,:)=correlation(Relax(7,:),n_sizerelax)

          Do i=1,n_sizerelax
             Write(8,'(i10,18(2x,E20.10E3))') int(Relax(1,i)),(Relax(j,i),j=2,18)
          enddo

          call SignatureFile(8,len_trim(fname),fname,len_trim('append'),'append')

          close(8)

    write(6,'(/,a,f8.4,2x,a,/)') 'Equilibrium files are written for T= ',kT/k_B,'Kelvin'

      endif


      END SUBROUTINE Relaxation
