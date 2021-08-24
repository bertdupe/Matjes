module m_llg_diag
use m_hamiltonian_collection, only: hamiltonian
use m_io_minimize,only : min_input
use m_constants,  only : hbar, mu_B
use m_io_files_utils
use m_io_utils
implicit none

private :: cart_to_sph, sph_to_cart, build_num_hess
public :: build_transmat 

interface diag_llg
	 module procedure build_transmat
end interface

contains


!build transition matrix
subroutine build_transmat(lat,io_simu,H)
    use m_derived_types, only : io_parameter,lattice
    type(io_parameter), intent(in)  :: io_simu
    type(lattice),intent(inout)     :: lat
    type(hamiltonian),intent(inout) :: H
    ! internal
    real(8), allocatable			:: M0(:,:)
	real(8),allocatable				:: Tr_theta(:,:), Tr_phi(:,:), Tr_thetaphi(:,:), Tr_phitheta(:,:), Tr(:,:)
	real(8),allocatable			 	:: Hess_theta(:,:), Hess_phi(:,:), Hess_thetaphi(:,:), Hess_phitheta(:,:)
	integer                     	:: N_cell,N_dim,N_mag
	integer							:: i,j,io
	real(8)							:: gp, hp, damping , gyro,mu_s  	
	logical 						:: print_tr,save_tr,save_tr_submat
	    
	! initialize
    N_cell=lat%Ncell
    N_dim=lat%M%dim_mode
    N_mag=(N_dim/3)*N_cell

    allocate(Tr_theta(N_mag,N_mag),Tr_phi(N_mag,N_mag),Tr_thetaphi(N_mag,N_mag),Tr_phitheta(N_mag,N_mag),Tr(2*N_mag,2*N_mag),source=0.0d0) 
    allocate(Hess_theta(N_mag,N_mag),Hess_phi(N_mag,N_mag),Hess_thetaphi(N_mag,N_mag),Hess_phitheta(N_mag,N_mag) ,source=0.0d0) 
	allocate(M0(2,N_mag),source=0.0d0) !sph
   
	io=open_file_read('input')
   	call get_parameter(io,'input','damping',damping)
   	gyro=1.0d0/(hbar)! if this is commented, frequencies in 1e15 Hz *1e-15) !=gamma/mu_s = mu_s/(hbar*mu_s)=1/hbar in 1/(eV.s), using seconds to get frequencies in Hz
   	gp=gyro/(1+damping**2) !precession term
   	hp=damping*gp !alignment term
 
   	print_tr=.false.
   	save_tr=.true.
   	save_tr_submat=.true.
   	   	
   	write(6,'(/,a,/)') 'Starting computation of the transition matrix of LLG. Warning: only for energy extrema.'
   	write(*,*) 'damping=', damping, ' gamma/mu_s= ',gyro,'1/(eV.fs)'
    write(*,*) 'gp=', gp, ' 1/(eV.fs), hp= ', hp, '1/(eV.fs)' 
    
    call cart_to_sph(lat,M0) !convert stable state to sph
   
	call build_num_hess(lat,io_simu,H,M0,Hess_theta,Hess_phi,Hess_thetaphi,Hess_phitheta) !build the hessian
	
    write(6,'(/,a,/)') 'Building the transition matrix...'
    
    Tr_theta=gp * Hess_phitheta - hp * Hess_theta
    
	do i=1,N_mag 
        do j=1,N_mag
    		Tr_phi(i,j)=- gp* sin(M0(1,j))/sin(M0(1,i)) * Hess_thetaphi(i,j) - hp* sin(M0(1,j))/sin(M0(1,i))  * Hess_phi(i,j)
    		Tr_thetaphi(i,j)= gp*sin(M0(1,j)) * Hess_phi(i,j) - hp *sin(M0(1,j))* Hess_thetaphi(i,j)
    		Tr_phitheta(i,j)= - gp/sin(M0(1,i)) * Hess_theta(i,j) - hp/sin(M0(1,i))* Hess_phitheta(i,j)
       	enddo
    enddo
    
    if(print_tr) then
			write(*,*)'Tr_theta='
		do j=1,N_mag  
		 	write(*,*)Tr_theta(j,:)
		 enddo
	   	 write(*,*)'Tr_phi='
	   	 do j=1,N_mag  
	   	 	write(*,*)Tr_phi(j,:)
	   	 enddo
	   	 write(*,*)'Tr_thetaphi='
	   	 do j=1,N_mag  
	   	 	write(*,*) Tr_thetaphi(j,:)
	   	 enddo
	   	 write(*,*)'Tr_phitheta='
	   	 do j=1,N_mag  
	   	 	write(*,*) Tr_phitheta(j,:)
		enddo
	end if
	
	if(save_tr) then
		!write to file
		open (1,file='trans_mat.dat')
		rewind 1
		do j=1,N_mag
	   		write(1,*) Tr_theta(j,:), Tr_thetaphi(j,:)
	   	enddo
		do j=1,N_mag
	   		write(1,*) Tr_phitheta(j,:), Tr_phi(j,:)
	   	enddo
	   	close(1)
	end if
	
	if(save_tr_submat) then
		!write to file
		open (1,file='Tr_theta.dat')
		open (2,file='Tr_thetaphi.dat')
		open (3,file='Tr_phitheta.dat')
		open (4,file='Tr_phi.dat')
		rewind 1
		rewind 2
		rewind 3
		rewind 4
		do j=1,N_mag
	   		write(1,*) Tr_theta(j,:)
	   		write(2,*) Tr_thetaphi(j,:)
	   		write(3,*) Tr_phitheta(j,:)
	   		write(4,*) Tr_phi(j,:)
	   	enddo
	   	close(1)
	   	close(2)
	   	close(3)
	   	close(4)
	end if
	
    
	write(6,'(/a,2x,E20.12E3/)') 'Done.'
	
end subroutine


! 					Build numerical Hessian matrix    
subroutine build_num_hess(lat,io_simu,H,M0,Hess_theta,Hess_phi,Hess_thetaphi,Hess_phitheta)
use m_derived_types, only : io_parameter,lattice
	use m_io_files_utils, only: open_file_write,close_file
    type(io_parameter), intent(in)  :: io_simu
    type(lattice),intent(inout)     :: lat
    type(hamiltonian),intent(inout) :: H
    real(8),intent(in) 				:: M0(:,:)
    real(8),intent(inout) 			:: Hess_theta(:,:), Hess_phi(:,:), Hess_thetaphi(:,:), Hess_phitheta(:,:)
    ! internal
    real(8) 						:: dphi, dtheta
    real(8), allocatable			:: Mp(:,:), Mm(:,:), Mpm(:,:), Mmp(:,:)!spherical coord
    real(8), allocatable	 		:: M0_cart(:,:),Mp_cart(:,:), Mm_cart(:,:), Mpm_cart(:,:), Mmp_cart(:,:) !cartesian coord
    real(8)							:: E0, Ep, Em, Epm, Emp !total energy
    integer                     	:: N_cell,N_dim,N_mag
    integer                     	:: i,j,io_file(2)
    real(8),pointer             	:: M3(:,:)
    logical :: print_hess,save_hess,save_hess_submat
    real(8)							:: epsi
    
    !initialization 
    
    write(6,'(/,a,/)') 'Building the numerical Hessian...'
    
    N_cell=lat%Ncell
    N_dim=lat%M%dim_mode
    N_mag=(N_dim/3)*N_cell
    
    allocate(Mp(2,N_mag),Mm(2,N_mag),Mpm(2,N_mag),Mmp(2,N_mag),source=0.0d0) !sph
    allocate(M0_cart(3,N_mag),Mp_cart(3,N_mag),Mm_cart(3,N_mag),Mpm_cart(3,N_mag),Mmp_cart(3,N_mag),source=0.0d0) !cart

    M3(1:3,1:N_mag)=>lat%M%all_modes
    M0_cart=M3

    Ep=0.0d0
    Em=0.0d0
    Epm=0.0d0
    Emp=0.0d0 
    
    dphi=0.001d0
    dtheta=0.001d0
    epsi=1.0d-10 !tolerance on accuracy
    
    print_hess=.false. 
    save_hess=.true. 
    save_hess_submat=.true. 
    
    E0=H%energy(lat)
    write(6,'(/a,2x,E20.12E3/)') 'Initial total energy density (eV/fu)',E0/real(N_cell,8)
    if(any(abs(sin(M0(1,:))).lt.1.0d-8)) stop 'This routine cannot be used if a spin lies at a pole.' !that's a problem
    
	Mp(:,:)=M0(:,:)
	Mm(:,:)=M0(:,:)
	Mmp(:,:)=M0(:,:)
	Mpm(:,:)=M0(:,:)

    do i=1,N_mag 
        do j=1,N_mag  
        	!--------- Hess_theta ---------!
			!small deviations
			Mp(1,i) = M0(1,i) + dtheta
			Mp(1,j) = M0(1,j) + dtheta
			
			Mm(1,i) = M0(1,i) - dtheta
			Mm(1,j) = M0(1,j) - dtheta
			
			Mpm(1,i)= M0(1,i) + dtheta
			Mpm(1,j)= M0(1,j) - dtheta
			
			Mmp(1,i)= M0(1,i) - dtheta
			Mmp(1,j)= M0(1,j) + dtheta
			
			!convert to cart
			call sph_to_cart(Mp,Mp_cart)
			call sph_to_cart(Mm,Mm_cart)
			call sph_to_cart(Mpm,Mpm_cart)	
			call sph_to_cart(Mmp,Mmp_cart)
	
			!compute energies
			M3=Mp_cart
			Ep=H%energy(lat)
			M3=Mm_cart
			Em=H%energy(lat)
			M3=Mmp_cart
			Emp=H%energy(lat)
			M3=Mpm_cart
			Epm=H%energy(lat)

			!Hess_theta entry
			if(i.ne.j) then !off-diagonal term
				Hess_theta(i,j)=(  Ep - Epm - Emp + Em ) / (4.0d0*dtheta*dtheta)
			else !diagonal term
				Hess_theta(i,j)= (  Ep - 2.0d0*E0 +  Em ) / (dtheta*dtheta) 
			endif	
			if(abs(Hess_theta(i,j)).le.epsi) Hess_theta(i,j)=0.0d0

			!revert to initial state
			Mp(:,:)=M0(:,:)
			Mm(:,:)=M0(:,:)
			Mmp(:,:)=M0(:,:)
			Mpm(:,:)=M0(:,:)		

        	!--------- Hess_phi ---------!
        	!small deviations
			Mp(2,i) = M0(2,i) + dphi
			Mp(2,j) = M0(2,j) + dphi
			
			Mm(2,i) = M0(2,i) - dphi
			Mm(2,j) = M0(2,j) - dphi
			
			Mpm(2,i)= M0(2,i) + dphi
			Mpm(2,j)= M0(2,j) - dphi
			
			Mmp(2,i)= M0(2,i) - dphi
			Mmp(2,j)= M0(2,j) + dphi
				
			!convert to cart
			call sph_to_cart(Mp,Mp_cart)
			call sph_to_cart(Mm,Mm_cart)
			call sph_to_cart(Mpm,Mpm_cart)	
			call sph_to_cart(Mmp,Mmp_cart)
				
			!compute energies
			M3=Mp_cart
			Ep=H%energy(lat)
			M3=Mm_cart
			Em=H%energy(lat)
			M3=Mmp_cart
			Emp=H%energy(lat)
			M3=Mpm_cart
			Epm=H%energy(lat)
			
			!Hess_phi entry
			if(i.ne.j) then !off-diagonal term
				Hess_phi(i,j)=(  Ep - Epm - Emp + Em ) /  ( 4.0d0 * dphi**2 *sin(M0(1,i))*sin(M0(1,j)) )
			else
				Hess_phi(i,j)= (  Ep - 2.0d0*E0 + Em ) / ( dphi**2 *sin(M0(1,i))*sin(M0(1,j)) )
			endif	
			if(abs(Hess_phi(i,j)).le.epsi) Hess_phi(i,j)=0.0d0
		
			!revert to initial state
			Mp(:,:)=M0(:,:)
			Mm(:,:)=M0(:,:)
			Mmp(:,:)=M0(:,:)
			Mpm(:,:)=M0(:,:)
			
			!--------- Hess_thetaphi ---------!
			!small deviations
			Mp(1,i) = M0(1,i) + dtheta
			Mp(2,j) = M0(2,j) + dphi
		
			Mm(1,i) = M0(1,i) - dtheta
			Mm(2,j) = M0(2,j) - dphi
			
			Mpm(1,i)= M0(1,i) + dtheta
			Mpm(2,j)= M0(2,j) - dphi
			
			Mmp(1,i)= M0(1,i) - dtheta
			Mmp(2,j)= M0(2,j) + dphi
				
			!convert to cart
			call sph_to_cart(Mp,Mp_cart)
			call sph_to_cart(Mm,Mm_cart)
			call sph_to_cart(Mpm,Mpm_cart)	
			call sph_to_cart(Mmp,Mmp_cart)
				
			!compute energies
			M3=Mp_cart
			Ep=H%energy(lat)
			M3=Mm_cart
			Em=H%energy(lat)
			M3=Mmp_cart
			Emp=H%energy(lat)
			M3=Mpm_cart
			Epm=H%energy(lat)
				
			!Hess_thetaphi entry
			Hess_thetaphi(i,j)=(  Ep - Epm - Emp + Em ) /  ( 4.0d0*dphi*dtheta*sin(M0(1,j)) )
			if(abs(Hess_thetaphi(i,j)).le.epsi) Hess_thetaphi(i,j)=0.0d0
			
			!revert to initial state
			Mp(:,:)=M0(:,:)
			Mm(:,:)=M0(:,:)
			Mmp(:,:)=M0(:,:)
			Mpm(:,:)=M0(:,:)
			
			!--------- Hess_phitheta ---------!
			!small deviations
			Mp(1,j) = M0(1,j) + dtheta
			Mp(2,i) = M0(2,i) + dphi
			
			Mm(1,j) = M0(1,j) - dtheta
			Mm(2,i) = M0(2,i) - dphi
			
			Mpm(1,j)= M0(1,j) + dtheta
			Mpm(2,i)= M0(2,i) - dphi
			
			Mmp(1,j)= M0(1,j) - dtheta
			Mmp(2,i)= M0(2,i) + dphi
				
			!convert to cart
			call sph_to_cart(Mp,Mp_cart)
			call sph_to_cart(Mm,Mm_cart)
			call sph_to_cart(Mpm,Mpm_cart)	
			call sph_to_cart(Mmp,Mmp_cart)
				
			!compute energies
			M3=Mp_cart
			Ep=H%energy(lat)
			M3=Mm_cart
			Em=H%energy(lat)
			M3=Mmp_cart
			Emp=H%energy(lat)
			M3=Mpm_cart
			Epm=H%energy(lat)
					
			!Hess_thetaphi entry
			Hess_phitheta(i,j)=(  Ep - Epm - Emp + Em ) /  ( 4.0*dphi*dtheta *sin(M0(1,i)) )
			if(abs(Hess_phitheta(i,j)).le.epsi)  Hess_phitheta(i,j)=0.0d0
			
			!revert to initial state
			Mp(:,:)=M0(:,:)
			Mm(:,:)=M0(:,:)
			Mmp(:,:)=M0(:,:)
			Mpm(:,:)=M0(:,:)
					
		enddo
	enddo
	M3=M0_cart
	
	write(6,'(/a,2x,E20.12E3/)') 'Done.'
	
	if(print_hess) then
		!print to console
		write(*,*)'Hess_theta='
		do j=1,N_mag  
		 	write(*,*)Hess_theta(j,:)
		 enddo
	   	 write(*,*)'Hess_phi='
	   	 do j=1,N_mag  
	   	 	write(*,*)Hess_phi(j,:)
	   	 enddo
	   	 
	   	 write(*,*)'Hess_thetaphi='
	   	 do j=1,N_mag  
	   	 	write(*,*) Hess_thetaphi(j,:)
	   	 enddo
	   	 write(*,*)'Hess_phitheta='
	   	 do j=1,N_mag  
	   	 	write(*,*) Hess_phitheta(j,:)
	   	 enddo	
   	 endif	
   	 
   	if(save_hess) then
		open (1,file='Hessian_mat.dat')
		rewind 1
		do i=1,N_mag
	   		write(1,*) Hess_theta(i,:), Hess_thetaphi(i,:)
	   	enddo
		do i=1,N_mag
	   		write(1,*) Hess_phitheta(i,:),Hess_phi(i,:)
	   	enddo
	   	close(1)
	end if
	
	if(save_hess_submat) then
		open (1,file='Hess_theta.dat')
		open (2,file='Hess_thetaphi.dat')
		open (3,file='Hess_phitheta.dat')
		open (4,file='Hess_phi.dat')
		rewind 1
		rewind 2
		rewind 3
		rewind 4
		do j=1,N_mag
	   		write(1,*) Hess_theta(j,:)
	   		write(2,*) Hess_thetaphi(j,:)
	   		write(3,*) Hess_phitheta(j,:)
	   		write(4,*) Hess_phi(j,:)
	   	enddo
	   	close(1)
	   	close(2)
	   	close(3)
	   	close(4)
	end if
	
end subroutine


!		Converts magnetization cart into spherical
subroutine cart_to_sph(lat,M_sph)
use m_derived_types, only : lattice
    type(lattice),intent(in)    		:: lat
    real(8)	, intent(out)			 	:: M_sph(:,:)
    !internal
    integer                     		:: N_cell,N_dim,N_mag
    integer                     		:: iomp
    real(8),pointer             		:: M3(:,:)
    
   
    N_cell=lat%Ncell
    N_dim=lat%M%dim_mode
    N_mag=(N_dim/3)*N_cell
    M3(1:3,1:N_mag)=>lat%M%all_modes
    
    !theta
    if(any(M3(3,:).gt.1.0d0)) stop 'some cart. components are >1, fix the cart_to_sph routine'
	M_sph(1,:)=acos(M3(3,:)) 
    
 	!phi 
	M_sph(2,:)=atan2(M3(2,:),M3(1,:)) 
		
end subroutine


!		Converts magnetization sph into cart
subroutine sph_to_cart(M_sph,M_cart)
    real(8),intent(in)		:: M_sph(:,:)
    real(8),intent(out) 	:: M_cart(:,:)
    integer     :: i

!changes to loop as cray compiler whines otherwise... internal errors...
    do i=1,size(M_sph,2)
        M_cart(1,i)=cos(M_sph(2,i)) * sin(M_sph(1,i))
        M_cart(2,i)=sin(M_sph(2,i)) * sin(M_sph(1,i))
        M_cart(3,i)=cos(M_sph(1,i))
    enddo

!   	M_cart(1,:)=cos(M_sph(2,:)) * sin(M_sph(1,:))
!   	M_cart(2,:)=sin(M_sph(2,:)) * sin(M_sph(1,:))
!   	M_cart(3,:)=cos(M_sph(1,:))
end subroutine

end module m_llg_diag
