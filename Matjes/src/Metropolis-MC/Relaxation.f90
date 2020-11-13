module m_relaxation
implicit none
contains
!
! ===============================================================
SUBROUTINE Relaxation(lat,io_MC,N_cell,E_total,E,Magnetization,qeulerp,qeulerm,kt,acc,rate,nb,cone,Hams,Q_neigh)
    use mtprng
    use m_Corre
    use m_constants, only : k_b
    use m_topocharge_all
    use m_store_relaxation
    use m_derived_types, only : point_shell_Operator,lattice
    use m_basic_types, only : vec_point
    use m_modes_variables, only : point_shell_mode
    use m_H_public
    use m_topo_commons
    use m_io_utils
    use m_io_files_utils
    use m_convert
    use m_MCstep
    use m_input_types,only: MC_input
    ! input
    type(lattice),intent(inout)     :: lat
    real(kind=8), intent(inout)     :: qeulerp,qeulerm,cone,acc,rate,E_total,magnetization(3),E(8)
    real(kind=8), intent(inout)     :: nb
    real(kind=8), intent(in)        :: kT
    type(MC_input),intent(in)       :: io_MC 
    integer, intent(in)             :: N_cell
    class(t_H), intent(in)          :: Hams(:)
    integer,intent(in)              :: Q_neigh(:,:) 
    ! a big table
    real(kind=8) :: Relax(18,io_MC%n_sizerelax),dumy(5)
    ! Slope Index
    integer :: i_relaxation,i_MC,io_relax
    ! dummy
    integer :: i,j,n_w_step
    character(len=50) :: fname
    
    write(6,'(/,a,/)') "starting relaxation"
    
    Relax=0.0d0
    n_w_step=io_MC%n_thousand/io_MC%n_sizerelax
    
    !     Monte Carlo steps for thermal equilibrium
    !     The first time it will take longer
    !     -----------------------------------------------------------------
    do i_relaxation=1,io_MC%n_thousand
    !         T_relax_1 is probably larger then T_relax, this is because
    !         the last step might take more time to relax than an
    !         the step to an unordered structure
                !Relaxation of the System
        Do i_MC=1,io_MC%T_relax*N_cell
            Call MCStep(lat,io_MC,N_cell,E_total,E,Magnetization,kt,acc,rate,nb,cone,Hams)
        enddo
        !In case T_relax set to zero at least one MCstep is done
        Call MCStep(lat,io_MC,N_cell,E_total,E,Magnetization,kt,acc,rate,nb,cone,Hams)
    
    ! calculate the topocharge
        dumy=get_charge(lat,Q_neigh)
        qeulerp=dumy(1)
        qeulerm=dumy(2)
    
    ! Write the Equilibrium files
        if (io_MC%print_relax) then
            if (mod(i_relaxation,n_w_step).eq.0) call store_relaxation(Relax,i_relaxation,dble(i_relaxation), &
                  &  sum(E)/dble(N_cell),E,dble(N_cell),kt,Magnetization,rate,cone,qeulerp,qeulerm)
        endif
    
    enddo   ! enddo over the relaxation loop
    !!!!!!!!!!!!!***************************!!!!!!
    !!!!!!!!!!!!!***************************!!!!!!
    
    write(6,'(/,a,f8.4,2x,a,/)') 'System is relaxed for T= ',kT/k_B,'Kelvin'
    
    !print the Equilibrium files
    if (io_MC%print_relax) then
    
    ! numbering of the files
        fname=convert('Equilibriumphi-',kT/k_B,'.dat')
        io_relax=open_file_write(fname)
    
        Write(io_relax,'(a)') "#   1:n_MC  2:Etot  3:T  4:M  5:rate  6:cone  7:Q  8:Exch  9:Zeeman  10:Ani &
       &         11:4S  12:DM  13:biq  14:dip  15:stoner  16:Chi(E)(t,T)  17:Chi(M)(t,T)  18:Chi(T)(t)"
    
        Relax(16,:)=correlation(Relax(2,:),io_MC%n_sizerelax)
        Relax(17,:)=correlation(Relax(4,:),io_MC%n_sizerelax)
        Relax(18,:)=correlation(Relax(7,:),io_MC%n_sizerelax)
    
        Do i=1,io_MC%n_sizerelax
             Write(io_relax,'(i10,18(2x,E20.10E3))') int(Relax(1,i)),(Relax(j,i),j=2,18)
        enddo
    
        call close_file(fname,io_relax)
    
        write(6,'(/,a,f8.4,2x,a,/)') 'Equilibrium files are written for T= ',kT/k_B,'Kelvin'
    
    endif

END SUBROUTINE Relaxation
end module
