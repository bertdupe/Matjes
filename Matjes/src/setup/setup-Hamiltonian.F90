module m_Hamiltonian
use m_derived_types
use m_io_utils
use m_io_files_utils
interface number_nonzero_coeff
  module procedure number_nonzero_coeff_1d
  module procedure number_nonzero_coeff_2d
end interface number_nonzero_coeff
private
public :: setup_Hamiltonians, get_Hamiltonians

contains

!!!!!!!!!!!!!!!!!!!!!!!
!  Routine that prepares the LOCAL Hamiltonian for a particular geometry
! The form of the Hamiltonian is
!
!  ( H_00  ...           )  (Site_0)
!  ( H_10  H_11          )  (Site_1)
!  (  .          .       )  (Site_2)
!  (  .             .    )  (...   )
!  (  .               .  )  (...   )
!
!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_Hamiltonians(Hamiltonian)
implicit none
type(Coeff_Ham), intent(out) :: Hamiltonian
! internal
integer :: io_param,neighbor_Jij,neighbor_DM,neighbor_ani
! maximum number of shell to consider
integer :: N_coeff_Exch,N_coeff_DMI,N_coeff_ani,N_coeff
type(site_Ham),allocatable,target,save :: shell_energy(:)
! exchange and DMI
real(kind=8), allocatable :: J_ij(:),D_ij(:),ani(:,:)
! slope
integer :: i,j

neighbor_Jij=0
neighbor_DM=0
neighbor_ani=0
io_param=open_file_read('input')

! count the neighbors for the exchange
neighbor_Jij=count_variables(io_param,'J_','input')
allocate(J_ij(neighbor_Jij))
J_ij=0.0d0

neighbor_DM=count_variables(io_param,'DM_','input')
allocate(D_ij(neighbor_DM))
D_ij=0.0d0

neighbor_ani=count_variables(io_param,'ani_','input')
allocate(ani(3,neighbor_ani))
ani=0.0d0

! get the exchange, DMI, anisotropy
call get_coeff(io_param,'input','J_',J_ij)
N_coeff_Exch=number_nonzero_coeff(J_ij,'exchange')
if (N_coeff_Exch.gt.0) Hamiltonian%i_exch=.true.

call get_coeff(io_param,'input','DM_',D_ij)
N_coeff_DMI=number_nonzero_coeff(D_ij,'DMI')
if (N_coeff_DMI.gt.0) Hamiltonian%i_DM=.true.

call get_parameter(io_param,'input','ani_',3,ani(:,1))
N_coeff_ani=number_nonzero_coeff(ani,'anisotropy')
if (N_coeff_ani.gt.0) Hamiltonian%i_ani=.true.

call close_file('input',io_param)

! total number of shell to consider
N_coeff=max(N_coeff_ani,N_coeff_DMI,N_coeff_Exch)+1

! allocate the variables
allocate(Hamiltonian%exchange(3,3,N_coeff_Exch))
allocate(Hamiltonian%ani(3,3,N_coeff_ani))
allocate(Hamiltonian%DMI(3,3,N_coeff_DMI))
allocate(Hamiltonian%total(3,3,N_coeff))
! initialize variables
Hamiltonian%exchange=0.0d0
Hamiltonian%ani=0.0d0
Hamiltonian%DMI=0.0d0
Hamiltonian%total=0.0d0

! set up the anisotropy
do i=1,N_coeff_ani
   do j=1,3
      Hamiltonian%ani(j,j,i)=ani(j,i)
   enddo
enddo

! set up the magnetic exchange interaction
do i=1,N_coeff_Exch
   do j=1,3
      Hamiltonian%exchange(j,j,i)=J_ij(i)
   enddo
enddo

! set up the DMI
do i=1,N_coeff_DMI
   Hamiltonian%DMI(2,1,i)=D_ij(i)
   Hamiltonian%DMI(3,1,i)=-D_ij(i)
   Hamiltonian%DMI(3,2,i)=D_ij(i)

   Hamiltonian%DMI(1,2,i)=-Hamiltonian%DMI(2,1,i)
   Hamiltonian%DMI(2,3,i)=-Hamiltonian%DMI(3,2,i)
   Hamiltonian%DMI(1,3,i)=-Hamiltonian%DMI(3,1,i)
enddo

! set up the total Hamiltonian per site
do i=1,N_coeff
   if (i.eq.1) then
      Hamiltonian%total(:,:,i)=Hamiltonian%ani(:,:,i)
   else
      if (i-1.le.N_coeff_Exch) Hamiltonian%total(:,:,i)=Hamiltonian%total(:,:,i)+Hamiltonian%exchange(:,:,i-1)
      if (i-1.le.N_coeff_DMI) Hamiltonian%total(:,:,i)=Hamiltonian%total(:,:,i)+Hamiltonian%DMI(:,:,i-1)
      if (i.le.N_coeff_ani) Hamiltonian%total(:,:,i)=Hamiltonian%total(:,:,i)+Hamiltonian%ani(:,:,i)
   endif
enddo

! allocate the target variable for the Hamiltonian and
! give values
allocate(shell_energy(N_coeff))
do i=1,N_coeff
  allocate(shell_energy(i)%H(3,3))
  shell_energy(i)%H(:,:)=Hamiltonian%total(:,:,i)
enddo

#ifdef CPP_DEBUG
do i=1,N_coeff
   write(6,*) shell_energy(i)%H(:,1)
   write(6,*) shell_energy(i)%H(:,2)
   write(6,*) shell_energy(i)%H(:,3)
enddo
#endif

end subroutine get_Hamiltonians

!
! function that determines how many coefficients are none 0 in the exchange or the DMI or whatever
!
integer function number_nonzero_coeff_2d(coeff,name)
implicit none
real(kind=8), intent(in) :: coeff(:,:)
character(len=*), intent(in) :: name
! internal
integer :: i,N(2)
real(kind=8) :: norm

N=shape(coeff(:,:))
number_nonzero_coeff_2d=0

! loop over the number of coefficients
do i=1,N(2)
   norm=sqrt(coeff(1,i)**2+coeff(2,i)**2+coeff(3,i)**2)
   if (norm.gt.1.0d-8) number_nonzero_coeff_2d=i
enddo

if (number_nonzero_coeff_2d.eq.0) then
   write(6,'(2a)') 'no non-zero coefficients found for ',name
endif

end function number_nonzero_coeff_2d

integer function number_nonzero_coeff_1d(coeff,name)
implicit none
real(kind=8), intent(in) :: coeff(:)
character(len=*), intent(in) :: name
! internal
integer :: i,N

N=size(coeff)
number_nonzero_coeff_1d=0

! loop over the number of coefficients
do i=1,N
   if (abs(coeff(i)).gt.1.0d-8) number_nonzero_coeff_1d=i
enddo

if (number_nonzero_coeff_1d.eq.0) then
   write(6,'(2a)') 'no non-zero coefficients found for ',name
endif

end function number_nonzero_coeff_1d

!!!!!!!!!!!!!!!!!!!!!!!
!  Routine that prepares the Hamiltonian for a particular geometry
!!!!!!!!!!!!!!!!!!!!!!!

      subroutine setup_Hamiltonians()
      use m_parameters
      use m_derived_types
      implicit none
! internal
      real(kind=8) :: energy,J1
      integer :: ic,il,i
      real(kind=8), save, target :: J_loc(3,3)
      real(kind=8), save, target :: zero(3,3)
      type(vec), save, target :: S(3)
      type(H_loc) :: voisin(3,3)
      type(vec_point) :: S2(3)
      logical :: test


      J_loc=0.0d0
      zero=0.0d0
      do i=1,3
         S(i)%w=(/0.0d0,0.0d0,1.0d0/)
         S2(i)%w => S(i)%w
      enddo

      J1=2.0d0

      J_loc(1,1)=1.0d0
      J_loc(2,2)=1.0d0
      J_loc(3,3)=1.0d0

      do ic=1,3
         do il=1,3
            nullify(voisin(il,ic)%H_loc)
         enddo
      enddo

      test = associated(voisin(1,2)%H_loc)
      if (.not.test) voisin(1,2)%H_loc => J_loc

      test = associated(voisin(2,1)%H_loc)
      if (.not.test) voisin(2,1)%H_loc => J_loc

      test = associated(voisin(2,3)%H_loc)
      if (.not.test) voisin(2,3)%H_loc => J_loc

      test = associated(voisin(3,2)%H_loc)
      if (.not.test) voisin(3,2)%H_loc => J_loc


!  make a routine of matrix multiply

      energy=0.0d0
      do ic=1,3
         do il=1,3
            test = associated(voisin(il,ic)%H_loc)
            if (.not.test) cycle
            energy = energy + get_E_local(S2(il)%w,voisin(il,ic)%H_loc,S(ic)%w,J1)
         enddo
      enddo

      end subroutine setup_Hamiltonians

      real(kind=8) function get_E_local(U,M,V,Jij)
      implicit none
      real(kind=8), intent(in) :: Jij
      real(kind=8), intent(in) :: V(3),U(3),M(3,3)
! internal
      real(kind=8) :: temp
      integer :: i,j

      get_E_local=0.0d0

      do i=1,3
         temp=0.0d0
         if (U(i).eq.0.0d0) cycle
         do j=1,3
            temp=temp+M(j,i)*V(j)
         enddo
         get_E_local=get_E_local+U(i)*temp*Jij
      enddo

      end function get_E_local

      end module m_Hamiltonian
