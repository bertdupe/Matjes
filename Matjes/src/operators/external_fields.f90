module m_external_fields
use m_derived_types, only : cell,vec_dim_n,lattice,vec_point
use m_operator_pointer_utils

type(vec_dim_N),private,protected,target,save :: EM_external

type(vec_point),public,protected,allocatable,save :: ext_field(:)

private
public :: get_EM_external_fields,associate_external_fields,initialize_external_fields

contains

!
! routine that establishes the dimension of the external electro-magnetic field
!
subroutine get_EM_external_fields(h_ext,E_ext,my_motif)
implicit none
real(kind=8), intent(in) :: h_ext(:),E_ext(:)
type(cell), intent(in) :: my_motif
! interenal variable
integer :: dim_H, dim_E,dim_EM
integer :: n_atom,n_mag,i,j
real(kind=8) :: dumy

dim_H=size(h_ext)
dim_E=size(E_ext)
dim_EM=0

!
! first find if there is an electric and a magnetic field
!

dumy=sqrt(sum(h_ext**2))
if (dumy.gt.1.0d-8) dim_EM=dim_EM+dim_H

dumy=sqrt(sum(E_ext**2))
if (dumy.gt.1.0d-8) dim_EM=dim_EM+dim_E

n_atom=size(my_motif%atomic)
n_mag=0

!
! first find the number of atom in the motif
!

do i=1,n_atom
   if (abs(my_motif%atomic(i)%moment).gt.1.0d-8) n_mag=n_mag+1
enddo

!
! Then the size of the EM_field is (dim(H)+dim(E))*n_mag
!

allocate(EM_external%w(dim_EM*n_mag))

EM_external%w=0.0d0

do i=1,n_mag
   if (abs(my_motif%atomic(i)%moment).lt.1.0d-8) cycle

   j=(i-1)*n_mag+1
   dumy=sqrt(sum(h_ext**2))
   if (dumy.gt.1.0d-8) EM_external%w(j:j+2)=h_ext*my_motif%atomic(i)%moment

   dumy=sqrt(sum(E_ext**2))
   if (dumy.gt.1.0d-8) EM_external%w(j+3:j+5)=E_ext

enddo

end subroutine get_EM_external_fields

!
! routine that associates the field pointer to the actual field value
!
subroutine associate_external_fields(tableNN)
implicit none
integer, intent(in) :: tableNN(:,:,:,:,:,:) !!tableNN(4,N_Nneigh,dim_lat(1),dim_lat(2),dim_lat(3),count(my_motif%i_mom)

call associate_pointer(ext_field,EM_external,tableNN)

end subroutine associate_external_fields

!
! routine initializes the EM_field depending on the number of atoms in the unit cell and so on
!

subroutine initialize_external_fields(my_lattice)
implicit none
type(lattice), intent(in) :: my_lattice
!internal variable
integer :: shape_lattice(4),size_point,i

shape_lattice=shape(my_lattice%l_modes)
size_point=product(shape_lattice)

allocate(ext_field(size_point))

do i=1,size_point
   nullify(ext_field(i)%w)
enddo

end subroutine initialize_external_fields

end module m_external_fields
