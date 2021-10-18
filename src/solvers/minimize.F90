module m_minimize
use m_H_public
use m_io_minimize,only : min_input
use m_derived_types, only : io_parameter,lattice
use m_hamiltonian_collection, only: hamiltonian
use mpi_basic
use m_work_ham_single, only:  work_ham_single
use, intrinsic :: iso_fortran_env, only : error_unit 
implicit none


private
public :: minimize, minimize_infdamp
public :: minimize_run, minimize_infdamp_run

contains



subroutine minimize(lat,io_simu,H,com)
    type(io_parameter), intent(in)  :: io_simu
    type(lattice), intent(inout)    :: lat
    type(hamiltonian),intent(inout) :: H
    type(mpi_type),intent(in)       :: com

    type(min_input)                 :: io_min

    if(com%ismas)then
        if(com%Np>1) write(error_unit,'(//A/A//)') "WARNING, USING MPI-PARALLELIZATION WHICH IS NOT IMPLEMENTED FOR MINIMIZE CALCULATION","this calculation will only run the the master thread" 
        Call io_min%read_file()
        Call minimize_run(lat,io_simu,io_min,H)
    endif
end subroutine


subroutine minimize_infdamp(lat,io_simu,H,com)
    type(io_parameter), intent(in)  :: io_simu
    type(lattice), intent(inout)    :: lat
    type(hamiltonian),intent(inout) :: H
    type(mpi_type),intent(in)       :: com

    type(min_input)                 :: io_min

    if(com%ismas)then
        if(com%Np>1) write(error_unit,'(//A/A//)') "WARNING, USING MPI-PARALLELIZATION WHICH IS NOT IMPLEMENTED FOR MINIMIZE_INFDAMP CALCULATION","this calculation will only run the the master thread" 
        Call io_min%read_file()
        Call minimize_infdamp_run(lat,io_simu,io_min,H)
    endif
end subroutine

!
!
! Minimization routine that it actually used.
! The interface is used to put the data into the good format
!
!
subroutine minimize_run(lat,io_simu,io_min,H)
    use m_derived_types, only : io_parameter,lattice
    use m_constants, only : pi
    use m_write_spin
    use m_createspinfile
    use m_solver, only : minimization
    use m_vector, only : norm_cross,norm, calculate_damping
    
    implicit none
    type(io_parameter), intent(in)  :: io_simu
    type(min_input), intent(in)     :: io_min
    type(lattice), intent(inout)    :: lat
    type(hamiltonian),intent(inout) :: H
    ! dummy variable
    real(8),allocatable, dimension(:,:)    :: velocity,predicator,force
    real(8),allocatable, dimension(:)      :: V_eff,F_temp
    real(8),allocatable,target             :: Feff(:)
    real(8),pointer     ::  Feff_vec(:,:)
    ! internal
    real(8)     :: dumy,force_norm,Energy,vmax,vtest,test_torque,max_torque
    ! the computation time
    integer(8)  :: i_min, gra_freq
    integer     :: gra_int
    logical :: gra_log
    integer :: iomp,dim_mode,N_cell

    gra_freq=int(io_simu%io_frequency,8)
    gra_log=io_simu%io_Xstruct
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    force_norm=0.0d0
    test_torque=0.0d0
    vmax=0.0d0
    vtest=0.0d0
    N_cell=lat%Ncell
    dim_mode=lat%M%dim_mode
    
    allocate(Feff(dim_mode*N_cell),source=0.0d0)
    Feff_vec(1:dim_mode,1:N_cell)=>Feff
    allocate(velocity(dim_mode,N_cell),predicator(dim_mode,N_cell),force(dim_mode,N_cell),source=0.0d0)
    allocate(V_eff(dim_mode),F_temp(dim_mode),source=0.0d0)
    
    Call H%get_eff_field(lat,Feff,1)

    do iomp=1,lat%Ncell
        force(:,iomp)=calculate_damping(lat%M%modes_v(:,iomp),Feff_vec(:,iomp))
        call minimization(lat%M%modes_v(:,iomp),force(:,iomp),predicator(:,iomp),io_min%dt**2,io_min%mass*2.0d0)
    enddo
    test_torque=norm_cross(predicator(:,iomp),force(:,iomp),1,3)

    lat%M%modes_v=predicator
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! end of initialization
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do i_min=1,io_min%N_minimization
        max_torque=0.0d0
        dumy=0.0d0
        force_norm=0.0d0
        vmax=0.0d0
        
        Call H%get_eff_field(lat,Feff,1)
        do iomp=1,N_cell
            F_temp=calculate_damping(lat%M%modes_v(:,iomp),Feff_vec(:,iomp))
            call minimization(velocity(:,iomp),(force(:,iomp)+F_temp)/2.0d0,V_eff,io_min%dt,io_min%mass)
            Feff_vec(:,iomp)=F_temp
            force(:,iomp)=Feff_vec(:,iomp)
            velocity(:,iomp)=V_eff
            dumy=dumy+dot_product(V_eff,Feff_vec(:,iomp))
            force_norm=force_norm+norm(Feff_vec(:,iomp))**2
        enddo
        
        if (abs(dumy).gt.1.0d-8) then
            do iomp=1,N_cell
                velocity(:,iomp)=dumy*force(:,iomp)/force_norm
            enddo
        else
            velocity=0.0d0
        endif
        
        do iomp=1,N_cell
            call minimization(lat%M%modes_v(:,iomp),velocity(:,iomp),force(:,iomp),predicator(:,iomp),io_min%dt,io_min%mass)
            test_torque=norm(force(:,iomp))
            vtest=norm(velocity(:,iomp))**2
            if (vtest.gt.vmax) vmax=vtest
            if (test_torque.gt.max_torque) max_torque=test_torque
        enddo
       
        lat%M%modes_v=predicator

        if (mod(i_min,io_min%Efreq).eq.0)then
            energy=H%energy(lat)
            write(6,'(/,a,2x,I20)') 'iteration',i_min
            write(6,'(a,2x,f14.11)') 'Energy of the system (eV/unit cell)',Energy/dble(N_cell)
            write(6,'(2(a,2x,f14.11,2x))') 'convergence criteria:',io_min%conv_torque,',Measured Torque:',max_torque
            write(6,'(a,2x,f14.11,/)') 'speed of displacements:',vmax
        endif
        
        if (gra_log.and.(mod(i_min-1,gra_freq).eq.0)) then
            gra_int=int((i_min-1)/int(gra_freq,8),4)
            call WriteSpinAndCorrFile(gra_int,lat%M%modes_v,'spin_minimization')
            call CreateSpinFile(gra_int,lat%M)
        endif
        
        if (io_min%conv_torque.gt.max_torque) then
            write(6,'(a)') 'minimization converged'
            exit
        endif
    enddo ! number of minimization steps
    nullify(Feff_vec)
end subroutine


subroutine minimize_infdamp_run(lat,io_simu,io_min,H)
    use m_derived_types, only : io_parameter,lattice
   ! use m_constants, only : pi
    use m_write_spin
    use m_createspinfile
    use m_vector, only : cross,norm
    type(io_parameter), intent(in)  :: io_simu
    type(min_input), intent(in)     :: io_min
    type(lattice),intent(inout)     :: lat
    type(hamiltonian),intent(inout) :: H
   ! type(work_ham_single),intent(inout) :: work ! idk what this is
    
    ! internal
    real(8)                     :: max_torque,test_torque,Edy
    integer(8)                  :: iter,gra_freq
    integer                     :: iomp
    logical                     :: gra_log
    integer                     :: gra_int
    integer                     :: N_cell,N_dim,N_mag
    !real(8),allocatable,target  :: F_eff(:)
    !real(8),allocatable         :: F_norm(:),torque(:)
    real(8),pointer             :: M3(:,:)!,F_eff3(:,:)
    logical                     :: conv
    real(8)					    :: dummy(3), Beff(3), torque(3), Beff_norm
    type(work_ham_single)       :: work !type containing work arrays for single energy evaluation
    
    write(6,'(/,a,/)') 'entering the infinite damping minimization routine'
    
    N_cell=lat%Ncell
    N_dim=lat%M%dim_mode
    N_mag=(N_dim/3)*N_cell
    
    !allocate(F_norm(N_mag),source=0.0d0)
    !allocate(F_eff(N_dim*N_cell),torque(N_dim*N_cell),source=0.0d0)
    !F_eff3(1:3,1:N_mag)=>F_eff
    M3(1:3,1:N_mag)=>lat%M%all_modes
    Call H%get_single_work(1,work)  !allocate work arrays for single energy evaluation (1 for magnetism)
    gra_log=io_simu%io_Xstruct
    gra_freq=int(io_simu%io_frequency,8)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Prepare the calculation of the energy and the effective field
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    Edy=H%energy(lat)
    write(6,'(/a,2x,E20.12E3/)') 'Initial total energy density (eV/fu)',Edy/real(N_cell,8)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !           Begin minimization
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !iter=0
    !max_torque=10.0d0
    do iter=1,io_min%N_minimization
        max_torque=0.0d0

        !Call H%get_eff_field(lat,F_eff,1)
       ! F_norm=norm2(F_eff3,1)
        !if(any(F_norm.lt.1.0d-8)) stop 'problem in the infinite damping minimization routine' !avoid divide by 0
       ! torque=cross(lat%M%all_modes,F_eff,1,N_dim*N_cell)
       ! test_torque=max(maxval(torque),-minval(torque))
      !  if ( abs(test_torque).gt.max_torque ) max_torque=test_torque
      
        !ITERATIVELY align the moments onto local normalized field. Recompute effective field at each site.
        do iomp=1,N_mag
        	!get local normalized field
         	Call H%get_eff_field_single(lat,iomp,Beff,work,1,dummy)
        	Beff_norm=norm(Beff)
        	if(Beff_norm.lt.1.0d-8) stop 'Beff=0, problem in the infinite damping minimization routine' !avoid dividing by 0
        	Beff=Beff/Beff_norm
        	!write(*,*) 'ok1'
        	!get local normalized torque and its largest component
        	torque=cross(M3(:,iomp),Beff,1,3)
        	test_torque=maxval(dabs(torque))
      	  	if (test_torque.gt.max_torque ) max_torque=test_torque
        	        	!write(*,*) 'ok2'
        	!align moment to field (unless m=0) 
        	if(norm(M3(:,iomp)).gt.1.0d-8)  M3(:,iomp)=Beff 
        	        	!write(*,*) 'ok3'	            
        enddo
        
!        iter=iter+1
        !print max_torque every io_min%Efreq iterations
        if (mod(iter,io_min%Efreq).eq.0) write(*,*) 'Max torque =',max_torque
    
        !write config to files
        if ((gra_log).and.(mod(iter,gra_freq).eq.0)) then
            gra_int=int((iter-1)/int(gra_freq,8),4)
            call WriteSpinAndCorrFile(gra_int,lat%M%modes_v,'spin_minimization')
            call CreateSpinFile(gra_int,lat%M)
            write(6,'(a,3x,I10)') 'wrote Spin configuration and povray file number',iter/gra_freq
        endif
        
        !test convergence
        conv=max_torque.lt.io_min%conv_torque
        if(conv)then
            write(*,*) 'Max_torque=',max_torque,' tolerance reached, minimization completed in ',iter,' iterations.'
            exit
        endif
    enddo
    if(.not.conv)then
        write(*,'(///A)') "WARNING, minimization routine did not reach minimium"
        write(*,*) 'Max_torque=            ',max_torque
        write(*,*) 'Convergence criterion= ',io_min%conv_torque
        write(*,'(///)')
    endif
    Edy=H%energy(lat)
    write(6,'(/a,2x,E20.12E3/)') 'Final total energy density (eV/fu)',Edy/real(N_cell,8)
    nullify(M3)
    !deallocate(F_eff,F_norm,torque)
end subroutine

end module m_minimize
