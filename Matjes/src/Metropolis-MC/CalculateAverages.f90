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
     qeulerp_av(i)=qeulerp_av(i)/Total_MC_Steps/pi(4.0d0)
     qeulerm_av(i)=qeulerm_av(i)/Total_MC_Steps/pi(4.0d0)
     chi_Q(1,i)=(Q_sq_sum_av(i)/Total_MC_Steps/pi(4.0d0)**2-(qeulerp_av(i)+qeulerm_av(i))**2)/kT(i)
     chi_Q(2,i)=(Qp_sq_sum_av(i)/Total_MC_Steps/pi(4.0d0)**2-qeulerp_av(i)**2)/kT(i)
     chi_Q(3,i)=(Qm_sq_sum_av(i)/Total_MC_Steps/pi(4.0d0)**2-qeulerm_av(i)**2)/kT(i)
     chi_Q(4,i)=((-Qm_sq_sum_av(i)/Total_MC_Steps/pi(4.0d0)**2-qeulerm_av(i)**2)* &
     &    (Qp_sq_sum_av(i)/Total_MC_Steps/pi(4.0d0)**2-qeulerp_av(i)**2))/kT(i)
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
qeulerp_av=qeulerp_av/Total_MC_Steps/pi(4.0d0)
qeulerm_av=qeulerm_av/Total_MC_Steps/pi(4.0d0)
chi_Q(1)=((qeulerp_av+qeulerm_av)**2-Q_sq_sum/Total_MC_Steps/pi(4.0d0)**2)/kT
chi_Q(2)=(qeulerp_av**2-Qp_sq_sum/Total_MC_Steps/pi(4.0d0)**2)/kT
chi_Q(3)=(qeulerm_av**2-Qm_sq_sum/Total_MC_Steps/pi(4.0d0)**2)/kT
chi_Q(4)=((-Qm_sq_sum/Total_MC_Steps/pi(4.0d0)**2-qeulerm_av**2)* &
     &   (Qm_sq_sum/Total_MC_Steps/pi(4.0d0)**2-qeulerp_av**2))/kT
vortex_av(:)=vortex_av(:)/Total_MC_Steps/3.0d0/sqrt(3.0d0)
if (Cor_log) chi_l(:)=total_MC_steps

END subroutine Calculate_thermo_serial

! ===============================================================
subroutine update_ave(sum_qp,sum_qm,Q_sq_sum,Qp_sq_sum,Qm_sq_sum,sum_vortex,vortex, &
     & E_sum,E_sq_sum,M_sum,M_sq_sum,E,Magnetization)
use m_topo_commons
Implicit none
real(kind=8), intent(in) :: E,Magnetization(:),vortex(:)
real(kind=8), intent(inout) ::sum_qm,sum_qp,sum_vortex(:),Q_sq_sum,Qp_sq_sum,Qm_sq_sum
real(kind=8), intent(inout) :: E_sum,E_sq_sum,M_sum(:),M_sq_sum(:)
! internal variables
real(kind=8) :: dumy(5),qeulerp,qeulerm

qeulerp=0.0d0
qeulerm=0.0d0
!     estimating the values M_av, E_av, S and C
!     and do so for any sublattice

E_sum=E_sum+E
E_sq_sum=E_sq_sum+E**2
! calculate the topocharge
dumy=get_charge()
qeulerp=dumy(1)
qeulerm=dumy(2)

sum_qp=sum_qp+qeulerp
sum_qm=sum_qm+qeulerm
Q_sq_sum=Q_sq_sum+(qeulerp+qeulerm)**2
Qp_sq_sum=Qp_sq_sum+qeulerp**2
Qm_sq_sum=Qm_sq_sum+qeulerm**2
M_sum=M_sum+Magnetization
M_sq_sum=M_sq_sum+Magnetization**2
sum_vortex=sum_vortex+dumy(3:5)
STOP 'update_ave'

END subroutine update_ave

! ===============================================================
subroutine initialize_ave(lat,Q_neigh,qeulerp,qeulerm,vortex,Magnetization)
!THIS RETURNS NO AVERAGES... IS THIS REALLY WHAT IS INTENDED TO BE CALCULATED?
    use m_topo_commons
    use m_derived_types
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
    qeulerm=dumy(2)
    vortex=dumy(3:5)
    nullify(M3)
END subroutine initialize_ave

end module m_average_MC
! ===============================================================
