module m_local_energy
use m_basic_types, only : vec_point
use m_derived_types, only : point_shell_Operator
use m_modes_variables, only : point_shell_mode
implicit none

interface local_energy
  module procedure  local_energy_pointer,local_energy_optimized
end interface

interface get_E_matrix
  module procedure  get_E_matrix_normal,get_E_matrix_T
end interface

! all vectors on one line (must be updated at each line)
real(kind=8), allocatable, dimension(:) :: all_vectors

! all E matrice on one line (must be done once since the table is always done in the same order)
real(kind=8), allocatable, dimension(:,:) :: all_E


#ifdef __sparse_mkl__
!sparse format full matrix
real(8),allocatable :: val(:)
integer,allocatable :: rowind(:),colind(:)
integer                         :: dimH
public :: set_H_sparse,energy_sparse,dimH
#endif

private
public :: local_energy,get_E_matrix,kill_E_matrix
#ifdef __direct_mult__
public :: set_large_H,energy_H
#endif

contains

!
! really too slow
!
subroutine local_energy_pointer(E_int,iomp,spin,dim_mode)
use m_energy_commons, only : energy
use m_dipole_energy
use m_dipolar_field, only : i_dip
use m_matrix, only : reduce
implicit none
! input
type(vec_point), intent(in) :: spin(:)
integer, intent(in) :: iomp,dim_mode
! ouput
real(kind=8), intent(out) :: E_int
! internal
integer :: i,N,j
real(kind=8) :: S_int(dim_mode)

N=size(energy%line(:,iomp))
E_int=0.0d0

do i=1,N

   j=energy%line(i,iomp)

   call reduce(energy%value(i,iomp),size(energy%value(i,iomp)%order_op),S_int,spin(iomp)%w,spin(j)%w,dim_mode)

   E_int=E_int+dot_product( spin(iomp)%w , S_int )

enddo

if (i_dip) E_int=E_int+get_dipole_E(iomp)

end subroutine local_energy_pointer


!
! much much faster
!

subroutine local_energy_optimized(E_int,iomp,spin)
use m_energy_commons, only : energy
use m_dipole_energy
use m_dipolar_field, only : i_dip
use m_matrix, only : reduce
#ifdef __EIGEN__
use m_eigen_interface, only: eigen_matmul_allE
#endif
implicit none
! input
type(vec_point), intent(in) :: spin(:)
integer, intent(in) :: iomp
! ouput
real(kind=8), intent(out) :: E_int
! internal
integer :: i,N,j,dim_mode

N=size(energy%line(:,iomp))
E_int=0.0d0
dim_mode=size(all_vectors)/N

!energy%line repr�sente toutes les lignes de H.
!Dans energy%line, la composante (i,iomp) correspond au voisin i
!du site iomp.
!Ici, i correspond aux valeurs non nulles
!j correspond � tous les voisins (j=1, N_voisins)
do i=1,N
   j=energy%line(i,iomp)
   all_vectors((i-1)*dim_mode+1:i*dim_mode)=spin(j)%w
enddo
#ifdef __EIGEN__
Call eigen_matmul_allE(size(spin(iomp)%w),spin(iomp)%w,size(all_vectors),all_vectors,E_int)
#else
E_int=dot_product( spin(iomp)%w , matmul( all_E , all_vectors ) )
#endif
if (i_dip) E_int=E_int+get_dipole_E(iomp)

end subroutine local_energy_optimized





subroutine get_E_matrix_normal(dim_mode)
use m_energy_commons, only : energy
use m_eigen_interface, only: eigen_set_all_E
implicit none
integer, intent(in) :: dim_mode
! internal
integer :: N,i,j

N=size(energy%line(:,1))

!
! check if the variables were already allocated
!

if (allocated(all_E).and.allocated(all_vectors)) return

allocate(all_vectors(dim_mode*N),all_E(dim_mode,dim_mode*N))
all_vectors=0.0d0
all_E=0.0d0

! the B_total is always read in the same direction so we can fill it only once
! one has to do a transpose here

do i=1,N
   all_E(:,(i-1)*dim_mode+1:i*dim_mode)=transpose(energy%value(i,1)%order_op(1)%Op_loc)
enddo

#ifdef __EIGEN__
Call eigen_set_all_E(size(all_E,1),size(all_E,2),all_E) !also do with get_E_matrix_T???
#endif
end subroutine get_E_matrix_normal

subroutine get_E_matrix_T(dim_mode,T)
use m_energy_commons, only : energy,translator
implicit none
integer, intent(in) :: dim_mode
real(kind=8), intent(in) :: T
! internal
integer :: N,i,j

N=size(energy%line(:,1))

! the B_total is always read in the same direction so we can fill it only once
! one has to do a transpose here

do i=1,N
   all_E(:,(i-1)*dim_mode+1:i*dim_mode)=transpose(energy%value(i,1)%order_op(1)%Op_loc)+translator(:,:,i)*T**2
enddo

end subroutine get_E_matrix_T



subroutine kill_E_matrix()
    implicit none
    deallocate(all_vectors,all_E)
    write(6,'(a)') 'Energy matrix deallocated'
end subroutine kill_E_matrix

#ifdef __sparse_mkl__
subroutine set_H_sparse(dim_mode)
    integer,intent(in)      ::  dim_mode

    Call set_matrix_sparse(dim_mode,val,rowind,colind,dimH)
end subroutine

subroutine energy_sparse(E,dimH)
    use m_lattice,only: modes
    real(8),intent(out) :: E
    integer,intent(in)  :: dimH
    integer             :: nnz
    real(8),pointer     :: mode(:)
    integer             :: N_site
    real(8)             :: tmp(dimH)
    external mkl_dcoogemv
   
    nnz=size(rowind)
    mode(1:dimH)=>modes
    Call mkl_dcoogemv('N',dimH,val,rowind,colind,nnz,mode,tmp)
    !Call mkl_dcoogemv('N',dimH,val,colind,rowind,nnz,mode,tmp)
    E=dot_product(tmp,mode)
end subroutine
#endif

subroutine set_matrix_sparse(dim_mode,val,rowind,colind,dimH)
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



#ifdef __direct_mult__
subroutine set_large_H(dim_mode)
    use m_eigen_interface, only: eigen_set_H
    integer,intent(in)  :: dim_mode
    real(8),allocatable :: val(:)
    integer,allocatable :: rowind(:),colind(:)
    integer             :: dimH

    Call set_matrix_sparse(dim_mode,val,rowind,colind,dimH)
    colind=colind-1
    rowind=rowind-1
    Call eigen_set_H(size(rowind),dimH,rowind,colind,val) 
end subroutine set_large_H


subroutine energy_H(E,dimmode)
    use m_lattice,only: modes
    use m_energy_commons, only : energy
    use m_eigen_interface, only: eigen_eval_H
    real(8),intent(out)     ::  E
    integer,intent(in)      ::  dimmode
    real(8),pointer         ::  mode(:)
    integer             :: N_site,dimH
    
    N_site=size(energy%line,2)
    dimH=N_site*dimmode
    mode(1:dimH)=>modes
    Call eigen_eval_H(dimH,mode,E)

end subroutine
#endif


end module m_local_energy
