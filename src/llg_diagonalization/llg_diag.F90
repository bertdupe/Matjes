module m_llg_diag
use m_H_public
use m_input_types,only : min_input
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
subroutine build_transmat(lat,io_simu,Hams)
    use m_derived_types, only : io_parameter,lattice
    type(io_parameter), intent(in)  :: io_simu
    type(lattice),intent(inout)     :: lat
    class(t_H), intent(in)          :: Hams(:)
    ! internal
    real(8), allocatable			:: M0(:,:)
	real(8),allocatable				:: Tr_theta(:,:), Tr_phi(:,:), Tr_thetaphi(:,:), Tr_phitheta(:,:), Tr(:,:)
	real(8),allocatable			 	:: Hess_theta(:,:), Hess_phi(:,:), Hess_thetaphi(:,:), Hess_phitheta(:,:)
	integer                     	:: N_cell,N_dim,N_mag
	integer							:: i,j,io
	real(8)							:: gp, hp, damping , gyro,mu_s  	
	logical 						:: print_tr,save_tr
	    
	! initialize
    N_cell=lat%Ncell
    N_dim=lat%M%dim_mode
    N_mag=(N_dim/3)*N_cell

    allocate(Tr_theta(N_mag,N_mag),Tr_phi(N_mag,N_mag),Tr_thetaphi(N_mag,N_mag),Tr_phitheta(N_mag,N_mag),Tr(2*N_mag,2*N_mag),source=0.0d0) 
    allocate(Hess_theta(N_mag,N_mag),Hess_phi(N_mag,N_mag),Hess_thetaphi(N_mag,N_mag),Hess_phitheta(N_mag,N_mag) ,source=0.0d0) 
	allocate(M0(2,N_mag),source=0.0d0) !sph
   
   	call get_parameter(io,'input','damping',damping)
   	gyro=1.0d0/hbar !1/(eV.s) 
   	!mu_s=2.7*mu_B !
   	gp=gyro/(1+damping**2) 
   	hp=damping*gp 
   	
   	print_tr=.false.
   	save_tr=.false.
   	   	
   	write(6,'(/,a,/)') 'Starting computation of the transition matrix of LLG. Warning: only for energy extrema.'
    
    call cart_to_sph(lat,M0) !convert stable state to sph
   
	call build_num_hess(lat,io_simu,Hams,M0,Hess_theta,Hess_phi,Hess_thetaphi,Hess_phitheta) !build the hessian
	
    write(6,'(/,a,/)') 'Building the transition matrix...'
    
    Tr_theta=gp * Hess_phitheta - hp * Hess_theta
    
	do i=1,N_mag 
        do j=1,N_mag
    		Tr_phi(i,j)=- gp* sin(M0(1,j))/sin(M0(1,i)) * Hess_thetaphi(i,j) - hp* sin(M0(1,j))/sin(M0(1,i))  * Hess_phi(i,j)
    		Tr_thetaphi(i,j)= gp*sin(M0(1,j)) * Hess_phi(i,j) - hp *sin(M0(1,j))* Hess_thetaphi(i,j);
    		Tr_phitheta(i,j)= - gp/sin(M0(1,i)) * Hess_theta(i,j) - hp/sin(M0(1,i))* Hess_phitheta(i,j);
       	enddo
    enddo
    
    if(print_tr) then
			write(*,*)'Tr_theta='
		do j=1,N_mag  
		 	write(*,*)Tr_theta(:,j)
		 enddo
	   	 write(*,*)'Tr_phi='
	   	 do j=1,N_mag  
	   	 	write(*,*)Tr_phi(:,j)
	   	 enddo
	   	 write(*,*)'Tr_thetaphi='
	   	 do j=1,N_mag  
	   	 	write(*,*) Tr_thetaphi(:,j)
	   	 enddo
	   	 write(*,*)'Tr_phitheta='
	   	 do j=1,N_mag  
	   	 	write(*,*) Tr_phitheta(:,j)
		enddo
	end if
	
	if(save_tr) then
		!write to file
		open (1,file='trans_mat.dat')
		rewind 1
		do j=1,N_mag
	   		write(1,*) Hess_theta(:,j), Hess_thetaphi(:,j)
	   	enddo
		do j=1,N_mag
	   		write(1,*) Hess_phitheta(:,j), Hess_phi(:,j)
	   	enddo
	   	close(1)
	end if
	
    
	write(6,'(/a,2x,E20.12E3/)') 'Done.'
	
end subroutine


! 					Build numerical Hessian matrix    
subroutine build_num_hess(lat,io_simu,Hams,M0,Hess_theta,Hess_phi,Hess_thetaphi,Hess_phitheta)
use m_derived_types, only : io_parameter,lattice
	use m_io_files_utils, only: open_file_write,close_file
    type(io_parameter), intent(in)  :: io_simu
   ! type(min_input), intent(in)     :: io_min
    type(lattice),intent(inout)     :: lat
    class(t_H), intent(in)          :: Hams(:)
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
    character(len=50)   :: file_name(2),form
    logical :: print_hess,save_hess
    
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
    
    print_hess=.true. 
    save_hess=.true. 
    
    E0=energy_all(Hams,lat)
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
			Ep=energy_all(Hams,lat)
			
			M3=Mm_cart
			Em=energy_all(Hams,lat)
			
			M3=Mmp_cart
			Emp=energy_all(Hams,lat)
			
			M3=Mpm_cart
			Epm=energy_all(Hams,lat)
			
			M3=M0

			!Hess_theta entry
			if(i.ne.j) then !off-diagonal term
				Hess_theta(i,j)=(  Ep - Epm - Emp + Em ) / (4.0d0*dtheta*dtheta)
				!write(*,*) 'i, j=',i,j,' Hess_theta(i,j)= ',Hess_theta(i,j),' where the num = ',(  Ep - Epm - Emp + Em ),' and the denom= ',(4.0d0*dtheta*dtheta)
				!write(*,*)'Ep= ', Ep, ' Em= ', Em, ' Epm= ',Epm !, ' Emp= ',Emp
				!write(*,*) 'for Mmp, theta_i,phi_i =', Mmp(1,i), Mmp(2,i), ' theta_j,phi_j= ',Mmp(1,j), Mmp(2,j)
				!write(*,*) 'Mmp_cart_i= ',Mmp_cart(:,i),'Mmp_cart_j=',Mmp_cart(:,j) 

				!write(*,*)' Emp= ',Emp  
			else !diagonal term
				Hess_theta(i,j)= (  Ep - 2.0d0*E0 + Em ) / (dtheta*dtheta)
				!write(*,*) 'i=j = ',i,j,' Hess_theta(i,j)= ',Hess_theta(i,j)
			endif	

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
			Ep=energy_all(Hams,lat)
			M3=Mm_cart
			Em=energy_all(Hams,lat)
			M3=Mmp_cart
			Emp=energy_all(Hams,lat)
			M3=Mpm_cart
			Epm=energy_all(Hams,lat)
			M3=M0
			
			!Hess_phi entry
			if(i.ne.j) then !off-diagonal term
				Hess_phi(i,j)=(  Ep - Epm - Emp + Em ) /  ( 4.0d0 * dphi**2 *sin(M0(1,i))*sin(M0(1,j)) )
			else
				Hess_phi(i,j)= (  Ep - 2.0d0*E0 + Em ) / ( dphi**2 *sin(M0(1,i))*sin(M0(1,i)) )
			endif	
			
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
			Ep=energy_all(Hams,lat)
			M3=Mm_cart
			Em=energy_all(Hams,lat)
			M3=Mmp_cart
			Emp=energy_all(Hams,lat)
			M3=Mpm_cart
			Epm=energy_all(Hams,lat)
			M3=M0
				
			!Hess_thetaphi entry
			Hess_thetaphi(i,j)=(  Ep - Epm - Emp + Em ) /  ( 4.0d0*dphi*dtheta*sin(M0(1,j)) )
			
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
			Ep=energy_all(Hams,lat)
			M3=Mm_cart
			Em=energy_all(Hams,lat)
			M3=Mmp_cart
			Emp=energy_all(Hams,lat)
			M3=Mpm_cart
			Epm=energy_all(Hams,lat)
			M3=M0	
		
			!Hess_thetaphi entry
			Hess_phitheta(i,j)=(  Ep - Epm - Emp + Em ) /  ( 4.0*dphi*dtheta *sin(M0(1,i)) )
			
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
	   	!write to file
	   	!file_name(1)='Hessian_mat.dat'
		!open(io_file(1), file=file_name(1))!open_file_write(file_name(1)) !crashes idk why
		open (1,file='Hessian_mat.dat')
		rewind 1
		do i=1,N_mag
	   		write(1,*) Hess_theta(i,:), Hess_thetaphi(i,:)
	   	enddo
		do i=1,N_mag
	   		write(1,*) Hess_phitheta(i,:),Hess_phi(i,:)
	   	enddo
	   	close(1)
	!	call close_file(file_name(1),io_file(1))
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
    
    M3(1:3,1:N_mag)=>lat%M%all_modes
    N_cell=lat%Ncell
    N_dim=lat%M%dim_mode
    N_mag=(N_dim/3)*N_cell
    
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
   	M_cart(1,:)=cos(M_sph(2,:)) * sin(M_sph(1,:))
   	M_cart(2,:)=sin(M_sph(2,:)) * sin(M_sph(1,:))
   	M_cart(3,:)=cos(M_sph(1,:))
end subroutine

end module m_llg_diag
