module m_internal_fields_commons
use m_derived_types, only : operator_real,Coeff_Ham,point_shell_Operator,vec_point,point_shell_mode
use m_operator_pointer_utils
! coefficients of the internal magnetic field
type(Coeff_Ham),public,protected,target,save :: B_int
! total effective field tensor
type(operator_real),public,protected,save :: B_total
! Exchange effective field tensor
type(operator_real),public,protected,save :: B_exchange
! DMI effective field tensor
type(operator_real),public,protected,save :: B_DMI
! anisotropy effective field tensor
type(operator_real),public,protected,save :: B_anisotropy
! Zeeman field tensor
type(operator_real),public,protected,save :: B_Zeeman


private
public :: associate_internal_Beff,associate_line_field

contains

!!!!!!!!!!!!!!!!!!!!!!!
!  Routine that prepares the LOCAL effective field based on the Hamiltonian in
!  energy_commons.f90
! The form of the effective field is the same as the Hamiltonian
!
!  ( B_00  ...           )  (Site_0)
!  ( B_10  B_11          )  (Site_1)
!  (  .          .       )  (Site_2)
!  (  .             .    )  (...   )
!  (  .               .  )  (...   )
!
!!!!!!!!!!!!!!!!!!!!!!!

subroutine associate_internal_Beff(my_lattice,tableNN,indexNN)
use m_derived_types, only : Coeff_Ham,lattice
use m_energy_commons, only : energy,Hamiltonian
use m_constants, only : identity
use m_get_position
implicit none
type(lattice), intent(in) :: my_lattice
integer, intent(in) :: tableNN(:,:,:,:,:,:) !!tableNN(4,N_Nneigh,dim_lat(1),dim_lat(2),dim_lat(3),count(my_motif%i_mom)
integer, intent(in) :: indexNN(:)
! internal
!type(operator_real),save :: energy
! slope of the sums
integer :: i,j,N_coeff_DMI,N_all_vois,N_all_shell,N_coeff_ani,N_coeff_Exch,l
integer :: Nl,Nc

Nl=energy%nline
Nc=energy%ncolumn

!Hamiltonian%total_shell(i)%atom(l)%H(2,1)
! allocate the variables
N_all_shell=size(Hamiltonian%total_shell)
N_all_vois=size(Hamiltonian%total_shell(1)%atom)
allocate(B_int%total_shell(N_all_shell))
do i=1,N_all_shell
   allocate(B_int%total_shell(i)%atom(N_all_vois))
   do j=1,N_all_vois
      allocate(B_int%total_shell(i)%atom(j)%H(3,3))
      B_int%total_shell(i)%atom(j)%H=0.0d0
   enddo
enddo

N_coeff_Exch=size(Hamiltonian%exchange,3)
allocate(B_int%exchange(3,3,N_coeff_Exch))
B_int%exchange=-2.0d0*Hamiltonian%exchange

N_coeff_ani=size(Hamiltonian%ani,3)
if (N_coeff_ani.eq.0) N_coeff_ani=1
allocate(B_int%ani(3,3,N_coeff_ani))
B_int%ani=-2.0d0*Hamiltonian%ani

N_coeff_DMI=size(Hamiltonian%DMI,3)
allocate(B_int%DMI(3,3,N_coeff_DMI))

do i=1,N_coeff_DMI
  B_int%DMI=Hamiltonian%DMI
enddo

! in case I want to put the magnetic field there
allocate(B_int%Zeeman(3,3))
B_int%Zeeman=-Hamiltonian%Zeeman
!B_int%total(:,:,1)=B_int%total(:,:,1)+B_int%Zeeman
!!!!!

if (N_coeff_ani.ne.0) B_int%total_shell(1)%atom(1)%H=B_int%total_shell(1)%atom(1)%H+B_int%ani(:,:,1)

do i=2,N_all_shell

   if (i-1.le.N_coeff_DMI) then

      if ( .not.all( abs(Hamiltonian%DMI(:,:,i-1)).lt.1.0d-8 ) ) then

         do l=1,indexNN(i-1)

            B_int%total_shell(i)%atom(l)%H=Hamiltonian%total_shell(i)%atom(l)%H
            B_int%total_shell(i)%atom(l)%H(1,1)=-2.0d0*Hamiltonian%total_shell(i)%atom(l)%H(1,1)
            B_int%total_shell(i)%atom(l)%H(2,2)=-2.0d0*Hamiltonian%total_shell(i)%atom(l)%H(2,2)
            B_int%total_shell(i)%atom(l)%H(3,3)=-2.0d0*Hamiltonian%total_shell(i)%atom(l)%H(3,3)

         enddo
      endif
   else

      if (i-1.le.N_coeff_Exch) B_int%total_shell(i)%atom(1)%H=B_int%total_shell(i)%atom(1)%H+B_int%exchange(:,:,i-1)

   endif

   if (i.le.N_coeff_ani) B_int%total_shell(i)%atom(1)%H=B_int%total_shell(i)%atom(1)%H+B_int%ani(:,:,i)
enddo

allocate(B_total%value(Nl,Nc))
B_total%nline=Nl
B_total%ncolumn=Nc

! first nullify all pointer
call dissociate(B_total,Nl,Nc)

call associate_pointer(B_total,B_int%total_shell,my_lattice,tableNN,indexNN)

#ifdef CPP_DEBUG
do i=1,N_all_shell
   write(*,*) 'shell', i

   do j=1,N_all_vois
      if (all( abs(Hamiltonian%total_shell(i)%atom(j)%H(:,:)).lt.1.0d-8 )) cycle

      write(6,*) B_int%total_shell(i)%atom(j)%H(1,:)
      write(6,*) B_int%total_shell(i)%atom(j)%H(2,:)
      write(6,*) B_int%total_shell(i)%atom(j)%H(3,:)
   enddo
enddo
#endif

end subroutine associate_internal_Beff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! optimize the calculation time
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine associate_line_field(N_cell,B_line,spin,mode_B_column)
implicit none
integer, intent(in) :: N_cell
type(vec_point),intent(in) :: spin(:)
type(point_shell_Operator), intent(inout) :: B_line(:)
type(point_shell_mode), intent(inout) :: mode_B_column(:)
!internal energy
integer :: i,k,j,n_atom_shell

do i=1,N_cell

   n_atom_shell=0
   do j=1,N_cell
      if (associated(B_total%value(j,i)%Op_loc)) n_atom_shell=n_atom_shell+1
   enddo

   allocate(B_line(i)%shell(n_atom_shell))
   allocate(mode_B_column(i)%shell(n_atom_shell))

   ! the k=1 case should be the case i=j
   ! it should be the first atom in the shell
   k=1
   do j=1,N_cell

      if (associated(B_total%value(j,i)%Op_loc)) then
         ! taking care of the on-site term
         if (i.eq.j) then
            B_line(i)%shell(1)%Op_loc=>B_total%value(j,i)%Op_loc
            mode_B_column(i)%shell(1)%w=>spin(j)%w
            cycle
         endif

         k=k+1
         B_line(i)%shell(k)%Op_loc=>B_total%value(j,i)%Op_loc
         mode_B_column(i)%shell(k)%w=>spin(j)%w
      endif

   enddo
enddo

end subroutine associate_line_field

end module m_internal_fields_commons
