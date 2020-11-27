module m_average_MC

interface CalculateAverages
  module procedure initialize_ave
  module procedure update_ave
end interface CalculateAverages

interface Calculate_thermo
  module procedure Calculate_thermo_serial
  module procedure Calculate_thermo_images
end interface Calculate_thermo


contains

! ===============================================================
subroutine Calculate_thermo_images(Cor_log,n_average, &
    &    N_cell,kT,E_sq_sum_av,E_sum_av,M_sq_sum_av, &
    &    C,chi_M,E_av,E_err_av,M_err_av,qeulerp_av,qeulerm_av,vortex_av,Q_sq_sum_av,Qp_sq_sum_av,Qm_sq_sum_av,chi_Q,chi_l, &
    &    M_sum_av,M_av,size_table)
use m_constants, only : pi
Implicit none
logical, intent(in) :: Cor_log
integer, intent(in) :: n_average,size_table
real(kind=8), intent(in) :: N_cell,kT(:),E_sq_sum_av(:),E_sum_av(:),M_sq_sum_av(:,:),M_sum_av(:,:),Q_sq_sum_av(:),Qp_sq_sum_av(:),Qm_sq_sum_av(:)
real(kind=8), intent(out) :: C(:),chi_M(:,:),E_av(:),chi_l(:,:),E_err_av(:),M_err_av(:,:),chi_Q(:,:),M_av(:,:)
real(kind=8), intent(inout) :: qeulerp_av(:),qeulerm_av(:),vortex_av(:,:)
! internal variables
real(kind=8) :: Total_MC_Steps
integer :: i

chi_l=0.0d0
E_err_av=0.0d0
M_err_av=0.0d0
Total_MC_Steps=dble(n_average)
chi_Q=0.0d0



do i=1,size_table
     C(i)=(E_sq_sum_av(i)/Total_MC_Steps-(E_sum_av(i)/(Total_MC_Steps))**2)/kT(i)**2/N_cell
     chi_M(:,i)=(M_sq_sum_av(:,i)/Total_MC_Steps-(M_sum_av(:,i)/Total_MC_Steps)**2)/kT(i)/N_cell
     if (n_average.gt.1) E_err_av(i)=sqrt(abs(E_sq_sum_av(i)-(E_sum_av(i))**2/Total_MC_Steps)/(Total_MC_Steps-1))/N_cell
     if (n_average.gt.1) M_err_av(:,i)=sqrt(abs(M_sq_sum_av(:,i)-M_sum_av(:,i)**2/Total_MC_Steps)/(Total_MC_Steps-1))/N_cell
     E_av(i)=E_sum_av(i)/Total_MC_Steps/N_cell
     M_av(:,i)=M_sum_av(:,i)/Total_MC_Steps/N_cell
     qeulerp_av(i)=qeulerp_av(i)/Total_MC_Steps/pi/4.0d0
     qeulerm_av(i)=qeulerm_av(i)/Total_MC_Steps/pi/4.0d0
     chi_Q(1,i)= (Q_sq_sum_av(i)/Total_MC_Steps/16.0d0/pi**2-(qeulerp_av(i)+qeulerm_av(i))**2)/kT(i) ! /!\ the sign of this changed compared to the serial subroutine?
     chi_Q(2,i)=(Qp_sq_sum_av(i)/Total_MC_Steps/16.0d0/pi**2-qeulerp_av(i)**2)/kT(i)
     chi_Q(3,i)=(Qm_sq_sum_av(i)/Total_MC_Steps/16.0d0/pi**2-qeulerm_av(i)**2)/kT(i)
     chi_Q(4,i)=((-Qm_sq_sum_av(i)/Total_MC_Steps/16.0d0/pi**2-qeulerm_av(i)**2)* &
     &    (Qp_sq_sum_av(i)/Total_MC_Steps/16.0d0/pi**2-qeulerp_av(i)**2))/kT(i)
     vortex_av(:,i)=vortex_av(:,i)/Total_MC_Steps/3.0d0/sqrt(3.0d0)
     if (Cor_log) chi_l(:,i)=total_MC_steps
enddo

END subroutine Calculate_thermo_images

! ===============================================================
subroutine Calculate_thermo_serial(Cor_log,n_average, &
    &    N_cell,kT,E_sq_sum_av,E_sum_av,M_sq_sum_av, &
    &    C,chi_M,E_av,E_err_av,M_err_av,qeulerp_av,qeulerm_av,vortex_av,Q_sq_sum,Qp_sq_sum,Qm_sq_sum, &
    &    chi_Q,chi_l, &
    &    M_sum_av,M_av)
    use m_constants, only : pi
    Implicit none
    logical, intent(in) :: Cor_log
    integer, intent(in) :: n_average
    real(kind=8), intent(in) :: N_cell,kT,E_sq_sum_av,E_sum_av,M_sq_sum_av(:),M_sum_av(:),Q_sq_sum,Qp_sq_sum,Qm_sq_sum
    real(kind=8), intent(out) :: C,chi_M(:),E_av,chi_l(:),E_err_av,M_err_av(:),chi_Q(:),M_av(:)
    real(kind=8), intent(inout) :: qeulerp_av,qeulerm_av,vortex_av(:)
    ! internal variables
    real(kind=8) :: Total_MC_Steps
    ! components of teh spins
    
    chi_l=0.0d0
    E_err_av=0.0d0
    M_err_av=0.0d0
    Total_MC_Steps=dble(n_average)
    chi_Q=0.0d0

    !      If compiled serial!
    C=(E_sq_sum_av/Total_MC_Steps-(E_sum_av/(Total_MC_Steps))**2)/kT**2/N_cell
    chi_M=(M_sq_sum_av(:)/Total_MC_Steps-(M_sum_av(:)/Total_MC_Steps)**2)/kT/N_cell
    if (n_average.gt.1) E_err_av=sqrt(abs(E_sq_sum_av-(E_sum_av)**2/Total_MC_Steps)/(Total_MC_Steps-1))/N_cell
    if (n_average.gt.1) M_err_av(:)=sqrt(abs(M_sq_sum_av(:)-M_sum_av(:)**2/Total_MC_Steps)/(Total_MC_Steps-1))/N_cell
    E_av=E_sum_av/Total_MC_Steps/N_cell
    M_av=M_sum_av(:)/Total_MC_Steps/N_cell
    qeulerp_av=qeulerp_av/Total_MC_Steps/pi/4.0d0
    qeulerm_av=qeulerm_av/Total_MC_Steps/pi/4.0d0

    chi_Q(1)=-((qeulerp_av-qeulerm_av)**2 - Q_sq_sum/Total_MC_Steps/16.0d0/pi**2)/kT/(qeulerp_av-qeulerm_av)
    chi_Q(2)=-(qeulerp_av**2-Qp_sq_sum/Total_MC_Steps/16.0d0/pi**2)/kT/qeulerp_av
    chi_Q(3)=-(qeulerm_av**2-Qm_sq_sum/Total_MC_Steps/16.0d0/pi**2)/kT/qeulerm_av
    chi_Q(4)=((Qm_sq_sum/Total_MC_Steps/16.0d0/pi**2-qeulerm_av**2)* &
         &   (Qp_sq_sum/Total_MC_Steps/16.0d0/pi**2-qeulerp_av**2))/kT/(qeulerm_av*qeulerp_av)

    vortex_av(:)=vortex_av(:)/Total_MC_Steps/3.0d0/sqrt(3.0d0)
    if (Cor_log) chi_l(:)=total_MC_steps

END subroutine Calculate_thermo_serial

! ===============================================================
subroutine update_ave(lat,Q_neigh,sum_qp,sum_qm,Q_sq_sum,Qp_sq_sum,Qm_sq_sum,sum_vortex,vortex, &
     & E_sum,E_sq_sum,M_sum,M_sq_sum,Re_MpMm_sum,Im_MpMm_sum,E,Magnetization)
    use m_topo_commons
    use m_derived_types,only: lattice
    use m_vector, only : cross
	use m_get_table_nn,only :get_table_nn
    Implicit none
    type(lattice),intent(in)  :: lat
    integer,intent(in)        :: Q_neigh(:,:)
    real(kind=8), intent(in) :: E,Magnetization(:),vortex(:)
    real(kind=8), intent(inout) ::sum_qm,sum_qp,sum_vortex(:),Q_sq_sum,Qp_sq_sum,Qm_sq_sum
    real(kind=8), intent(inout) :: E_sum,E_sq_sum,M_sum(:),M_sq_sum(:)
    real(kind=8), intent(inout) ::  Re_MpMm_sum(:), Im_MpMm_sum(:)
    ! internal variables
    real(kind=8) :: dumy(5),qeulerp,qeulerm, temp(3)
    real(kind=8) :: Mi(3), Mj(3), Re_MpMm, Im_MpMm
    integer      :: N_neighbours,N_cell
    integer, allocatable    :: indexNN(:),tableNN(:,:,:,:,:,:)
    integer      :: i_cell,i_nei
    real(8),pointer :: M3(:,:)

    !     estimating the values M_av, E_av, S and C
    !     and do so for any sublattice
    
    E_sum=E_sum+E
    E_sq_sum=E_sq_sum+E**2

    ! calculate the topocharge
    dumy=get_charge(lat,Q_neigh)
    qeulerp=dumy(1)
    qeulerm=-dumy(2)
    
    sum_qp=sum_qp+qeulerp
    sum_qm=sum_qm+qeulerm
    Q_sq_sum=Q_sq_sum+(qeulerp-qeulerm)**2
    Qp_sq_sum=Qp_sq_sum+qeulerp**2
    Qm_sq_sum=Qm_sq_sum+qeulerm**2
    M_sum=M_sum+Magnetization
    M_sq_sum=M_sq_sum+Magnetization**2
	sum_vortex=sum_vortex+dumy(3:5)


    ! ------- for <M+M-> ------- !
    Re_MpMm=0.0d0
    Im_MpMm=0.0d0
   
    M3(1:3,1:lat%nmag*product(lat%dim_lat))=>lat%M%all_modes

    !get neighbours
    N_neighbours = 1
    Call get_table_nn(lat,N_neighbours,indexNN,tableNN)
    N_cell = lat%Ncell
    
	!write(*,*)"size of indexNN= ", size(indexNN), "shape of tableNN = ", shape(tableNN) 
	!write(*,*)"content of indexNN = ", indexNN
	!write(*,*)"content of tableNN = "
    !write(*,*) tableNN(:,:,:,:,:,:)
	!stop

	do i_cell=1,N_cell 
		Mi=M3(:,i_cell)
		do i_nei=1,N_neighbours
			nei=tableNN(1,2,3,4,5,6)
			Mj=M3(:,i_nei)
			Re_MpMm = Re_MpMm + dot_product(Mi,Mj) - Mi(3)*Mj(3) !real part M+M-, sum over all pairs
			temp = cross(Mi,Mj)
			Im_MpMm= Im_MpMm + temp(3) !imaginary part of M+M-, sum over all pairs
		enddo
	enddo

	Re_MpMm_sum=Re_MpMm_sum+Re_MpMm
	Im_MpMm_sum=Im_MpMm_sum+Im_MpMm

END subroutine update_ave

! ===============================================================
subroutine initialize_ave(lat,Q_neigh,qeulerp,qeulerm,vortex,Magnetization)
!THIS RETURNS NO AVERAGES... IS THIS REALLY WHAT IS INTENDED TO BE CALCULATED?
    use m_topo_commons
    use m_derived_types,only: lattice
    Implicit none
    type(lattice),intent(in)  :: lat
    integer,intent(in)        :: Q_neigh(:,:)
    real(kind=8), intent(out) :: vortex(3),qeulerp,qeulerm,Magnetization(3)
    ! dumy
    real(kind=8) :: dumy(5)
    real(8),pointer :: M3(:,:)

    M3(1:3,1:lat%nmag*product(lat%dim_lat))=>lat%M%all_modes
    Magnetization=sum(M3,2) !sum over magnetization of all magnetic atoms without magnetic moment, probably not what is really wanted
    dumy=get_charge(lat,Q_neigh)
    qeulerp=dumy(1)
    qeulerm=-dumy(2)
    vortex=dumy(3:5)

    nullify(M3)
END subroutine initialize_ave

end module m_average_MC
! ===============================================================

!put this in different file later
subroutine get_neighbours(lat,flat_nei)
    use m_derived_types,only: lattice
	use m_get_table_nn,only :get_table_nn
    Implicit none
    type(lattice),intent(in)  :: lat
    integer,intent(inout)        :: flat_nei(:,:)
    ! internal variables
    integer      :: N_neighbours,N_cell
    integer, allocatable    :: indexNN(:),tableNN(:,:,:,:,:,:)

 
    N_neighbours = 1 !first neighbours
    Call get_table_nn(lat,N_neighbours,indexNN,tableNN)

	!convert to a table with linear indices (Ncell, indexNN)
	do i_z=1,shape_tableNN(5)
    	do i_y=1,shape_tableNN(4)
     	   do i_x=1,shape_tableNN(3)
				!if(tableNN(5,i_vois,i_x,i_y,i_z,1)/=1) cycle
                        !N_line=N_line+1
				tableNN(1:3,,i_x,i_y,i_z,1)/=1
                i_flat=[i_x,i_y,i_z]
                        !connect(1,N_line)=lat%index_m_1(ilat_1)
                        ilat_2=tableNN(1:3,i_vois+offset,i_x,i_y,i_z,1)
                        !connect(2,N_line)=lat%index_m_1(ilat_2)
	ii1=lat%index_m_1(ilat_1)


end subroutine get_neighbours
