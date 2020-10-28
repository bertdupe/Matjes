#ifdef CPP_MKL_SPBLAS
module mkl_spblas_util
use MKL_SPBLAS, only: mkl_sparse_z_export_csr,SPARSE_MATRIX_T
USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_DOUBLE_COMPLEX,C_PTR,C_F_POINTER,C_INT
private
public :: unpack_csr_to_coo,unpack_csr

contains


   subroutine unpack_csr(dimH,H,nnz,ia,ja,acsr)
        type(SPARSE_MATRIX_T),intent(in)   :: H
        !fortran csr output
        integer,intent(in)                              :: dimH
        complex(C_DOUBLE_COMPLEX),pointer,intent(out)   :: acsr(:)
        integer(C_INT),pointer,intent(out)              :: ia(:),ja(:)
        integer,intent(out)                             :: nnz


        !export csr parameters
        integer(C_INT)                     :: nrows,ncols,indexing
        integer(c_int),pointer             :: tmp(:)
        !integer                           :: nnz
        integer                            :: stat
        type(C_PTR)                        :: values
        type(C_PTR)                        :: rows_start,rows_end,col_indx

        !get csr matrix out of SPARSE_MATRIX_T
        STOP "UPDATE unpack_csr to also work with potentially different mkl interface, so far comment out at TB not immediately necessary"
        !stat=mkl_sparse_z_export_csr ( H , indexing , nrows , ncols , rows_start , rows_end , col_indx , values ) 
        if(stat/=0) STOP 'sparse_z_export_csr in uppack_csr failed'
        CAll C_F_POINTER(rows_end, tmp,[nrows] )
        nnz=tmp(nrows)-indexing
        CAll C_F_POINTER(values, acsr,[nnz]) 
        CAll C_F_POINTER(col_indx, ja,[nnz]) 
        CAll C_F_POINTER(rows_start, ia,[dimH+1])
   end subroutine


   subroutine unpack_csr_to_coo(dimH,H,nnz,rowind,colind,acoo)
        type(SPARSE_MATRIX_T),intent(in)   :: H
        integer,intent(in)                 :: dimH
        !convert to coordinate representation
        complex(8),allocatable,intent(out) :: acoo(:)
        integer,allocatable,intent(out)    :: rowind(:),colind(:)
        integer,intent(out)                :: nnz
        !fortran csr output 
        integer                            :: job(6)
        complex(C_DOUBLE_COMPLEX),pointer  :: acsr(:)
        integer(C_INT),pointer             :: ia(:),ja(:)

        integer                            :: info


        Call unpack_csr(dimH,H,nnz,ia,ja,acsr)

        !create coo sparse matrix and add negative branch in (2,2) quadrant
        job=[0,1,1,0,nnz,3]
        allocate(acoo(nnz),source=cmplx(0.0d0,0.0d0,8))
        allocate(rowind(nnz),source=0)
        allocate(colind(nnz),source=0)
        Call mkl_zcsrcoo ( job , dimH , acsr , ja , ia , nnz , acoo , rowind , colind , info )
        if(info/=0) STOP 'info/=0 in unpack_csr_to_coo'

   end subroutine
end module
#endif
