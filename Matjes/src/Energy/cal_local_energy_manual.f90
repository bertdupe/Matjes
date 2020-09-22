module m_local_energy_manual
use m_basic_types, only : vec_point
use m_derived_types, only : point_shell_Operator
use m_modes_variables, only : point_shell_mode
implicit none

! all vectors on one line (must be updated at each line)
real(kind=8), allocatable, dimension(:) :: all_vectors

! all E matrice on one line (must be done once since the table is always done in the same order)
real(kind=8), allocatable, dimension(:,:) :: all_E

private
public :: local_energy_manual,sum_energy_manual,set_E_matrix_manual,kill_E_matrix_manual,get_matrix_sparse_manual

contains


subroutine sum_energy_manual(E_int,lat)
use m_derived_types, only: lattice
implicit none
type(lattice), intent(in) :: lat
real(8), intent(out) :: E_int
! internal
real(8) :: E_i
integer :: i,N,j

E_int=0.0d0
do i=1,product(lat%dim_lat)
    call local_energy_manual(E_i,i,lat)
    E_int=E_int+E_i
enddo

end subroutine 


subroutine local_energy_manual(E_int,iomp,lat)
use m_energy_commons, only : energy
use m_dipole_energy
use m_dipolar_field, only : i_dip
use m_derived_types, only: lattice
implicit none
! input
type(lattice), intent(in) :: lat
integer, intent(in) :: iomp
! output
real(kind=8), intent(out) :: E_int
! internal
integer :: i,N,j

N=size(energy%line(:,iomp))
E_int=0.0d0

!energy%line represente toutes les lignes de H.
!Dans energy%line, la composante (i,iomp) correspond au voisin i
!du site iomp.
!Ici, i correspond aux valeurs non nulles
!j correspond a tous les voisins (j=1, N_voisins)
do i=1,N
   j=energy%line(i,iomp)
   all_vectors((i-1)*lat%dim_mode+1:i*lat%dim_mode)=lat%ordpar%all_l_modes(j)%w
enddo

E_int=dot_product(lat%ordpar%all_l_modes(iomp)%w , matmul( all_E , all_vectors ) )

if (i_dip) E_int=E_int+get_dipole_E(iomp)

end subroutine 


subroutine set_E_matrix_manual(dim_mode)
use m_energy_commons, only : energy
use m_eigen_interface, only: eigen_set_all_E
implicit none
integer, intent(in) :: dim_mode
! internal
integer :: N,i,j

!
! check if the variables were already allocated
!

if (allocated(all_E).and.allocated(all_vectors)) return

N=size(energy%line(:,1))
allocate(all_vectors(dim_mode*N),all_E(dim_mode,dim_mode*N))
all_vectors=0.0d0
all_E=0.0d0

! the B_total is always read in the same direction so we can fill it only once
! one has to do a transpose here

do i=1,N
   all_E(:,(i-1)*dim_mode+1:i*dim_mode)=transpose(energy%value(i,1)%order_op(1)%Op_loc)
enddo

end subroutine 


subroutine kill_E_matrix_manual()
    implicit none
    deallocate(all_vectors,all_E)
    write(6,'(a)') 'Energy matrix deallocated'
end subroutine 


subroutine get_matrix_sparse_manual(dim_mode,val,rowind,colind,dimH)
    use m_energy_commons, only : energy
    integer,intent(in)      ::  dim_mode
    real(8),intent(out),allocatable :: val(:)
    integer,intent(out),allocatable :: rowind(:),colind(:)
    integer,intent(out)     :: dimH
    
    integer             :: nnz
    integer             :: i,j,l,ii
    integer             :: i1,i2
    integer             :: N_neigh
    integer             :: N_persite,N_site
    
    N_neigh=size(energy%line,1)
    N_site=size(energy%line,2)
    dimH=N_site*dim_mode
    N_persite=0
    do i=1,N_neigh
        N_persite=N_persite+count(energy%value(i,1)%order_op(1)%Op_loc /= 0.0d0)
    enddo
    nnz=N_persite*N_site
    allocate(val(nnz),source=0.0d0)
    allocate(colind(nnz),source=0)
    allocate(rowind(nnz),source=0)
    ii=1
    do i=1,N_site
        do l=1,N_neigh
            j=energy%line(l,i)
            do i2=1,dim_mode
                do i1=1,dim_mode
                    if(energy%value(l,i)%order_op(1)%Op_loc(i1,i2)/= 0.0d0)then
                        colind(ii)=(j-1)*dim_mode+i1
                        rowind(ii)=(i-1)*dim_mode+i2
                        val(ii)=energy%value(l,i)%order_op(1)%Op_loc(i1,i2)
                        ii=ii+1
                    endif
                enddo
            enddo
        enddo
    enddo
    
end subroutine
end module 
