#ifdef CPP_MKL
module mkl_spblas_util
use MKL_SPBLAS, only: mkl_sparse_z_export_csr,SPARSE_MATRIX_T
USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_DOUBLE_COMPLEX,C_PTR,C_F_POINTER,C_INT, C_DOUBLE
private
public :: unpack_csr_to_coo,unpack_csr

interface unpack_csr
    module procedure unpack_csr_complex
    module procedure unpack_csr_real
end interface

interface
    !manually write interface again since these interfaces seem to differ between different mkl versions
    function mkl_sparse_d_export_csr_cf(source,indexing,rows,cols,rows_start,rows_end,col_indx,values) bind(c, name='MKL_SPARSE_D_EXPORT_CSR')
        use iso_c_binding
        import sparse_matrix_t
        type(sparse_matrix_t) , intent(in) :: source
        integer(c_int), intent(inout) :: indexing
        integer(c_int) , intent(inout) :: rows
        integer(c_int) , intent(inout) :: cols
        type(c_ptr) , intent(inout) :: rows_start
        type(c_ptr) , intent(inout) :: rows_end
        type(c_ptr) , intent(inout) :: col_indx
        type(c_ptr) , intent(inout) :: values
        integer(c_int) mkl_sparse_d_export_csr_cf
    end function

    function mkl_sparse_z_export_csr_cf(source,indexing,rows,cols,rows_start,rows_end,col_indx,values) bind(c, name='MKL_SPARSE_Z_EXPORT_CSR')
        use iso_c_binding
        import sparse_matrix_t
        type(sparse_matrix_t) , intent(in) :: source
        integer(c_int), intent(inout) :: indexing
        integer(c_int) , intent(inout) :: rows
        integer(c_int) , intent(inout) :: cols
        type(c_ptr) , intent(inout) :: rows_start
        type(c_ptr) , intent(inout) :: rows_end
        type(c_ptr) , intent(inout) :: col_indx
        type(c_ptr) , intent(inout) :: values
        integer(c_int) mkl_sparse_z_export_csr_cf
    end function
end interface


contains
   subroutine unpack_csr_real(H,nnz,ia,ja,acsr,dim_mat)
        type(SPARSE_MATRIX_T),intent(in)   :: H
        !fortran csr output
        real(C_DOUBLE),pointer,intent(out)              :: acsr(:)
        integer(C_INT),pointer,intent(out)              :: ia(:),ja(:)
        integer,intent(out)                             :: nnz
        integer,intent(out),optional                    :: dim_mat(2)

        !export csr parameters
        integer(C_INT)                     :: nrows,ncols,indexing
        integer(c_int),pointer             :: tmp(:)
        !integer                           :: nnz
        integer                            :: stat
        type(C_PTR)                        :: values
        type(C_PTR)                        :: rows_start,rows_end,col_indx

        !get csr matrix out of SPARSE_MATRIX_T
        stat=mkl_sparse_d_export_csr_cf ( H , indexing , nrows , ncols , rows_start , rows_end , col_indx , values ) 
        if(stat/=0) STOP 'sparse_z_export_csr in unpack_csr failed'
        CAll C_F_POINTER(rows_end, tmp,[nrows] )
        nnz=tmp(nrows)-indexing
        CAll C_F_POINTER(values, acsr,[nnz]) 
        CAll C_F_POINTER(col_indx, ja,[nnz]) 
        CAll C_F_POINTER(rows_start, ia,[nrows+1])
        if(present(dim_mat)) dim_mat=[nrows,ncols]
   end subroutine

   subroutine unpack_csr_complex(H,nnz,ia,ja,acsr,dim_mat)
        type(SPARSE_MATRIX_T),intent(in)   :: H
        !fortran csr output
        complex(C_DOUBLE_COMPLEX),pointer,intent(out)   :: acsr(:)
        integer(C_INT),pointer,intent(out)              :: ia(:),ja(:)
        integer,intent(out)                             :: nnz
        integer,intent(out),optional                    :: dim_mat(2)

        !export csr parameters
        integer(C_INT)                     :: nrows,ncols,indexing
        integer(c_int),pointer             :: tmp(:)
        integer                            :: stat
        type(C_PTR)                        :: values
        type(C_PTR)                        :: rows_start,rows_end,col_indx

        !get csr matrix out of SPARSE_MATRIX_T
        stat=mkl_sparse_z_export_csr_cf ( H , indexing , nrows , ncols , rows_start , rows_end , col_indx , values ) 
        if(stat/=0) STOP 'sparse_z_export_csr in uppack_csr failed'
        CAll C_F_POINTER(rows_end, tmp,[nrows] )
        nnz=tmp(nrows)-indexing
        CAll C_F_POINTER(values, acsr,[nnz]) 
        CAll C_F_POINTER(col_indx, ja,[nnz]) 
        CAll C_F_POINTER(rows_start, ia,[nrows+1])
        if(present(dim_mat)) dim_mat=[nrows,ncols]
   end subroutine


    subroutine unpack_csr_to_coo(H,nnz,rowind,colind,acoo)
        !obsolete?    
        type(SPARSE_MATRIX_T),intent(in)   :: H
        !convert to coordinate representation
        complex(8),allocatable,intent(out) :: acoo(:)
        integer,allocatable,intent(out)    :: rowind(:),colind(:)
        integer,intent(out)                :: nnz
        !fortran csr output 
        integer                            :: job(6)
        complex(C_DOUBLE_COMPLEX),pointer  :: acsr(:)
        integer(C_INT),pointer             :: ia(:),ja(:)
        
        integer                            :: dim_mat(2)
        integer                            :: info

        Call unpack_csr(H,nnz,ia,ja,acsr,dim_mat)
        if(dim_mat(1)/=dim_mat(2)) ERROR STOP "Cannot unpack csr to coo for asymmetric shape"

        !create coo sparse matrix and add negative branch in (2,2) quadrant
        job=[0,1,1,0,nnz,3]
        allocate(acoo(nnz),source=cmplx(0.0d0,0.0d0,8))
        allocate(rowind(nnz),source=0)
        allocate(colind(nnz),source=0)
        Call mkl_zcsrcoo ( job , dim_mat , acsr , ja , ia , nnz , acoo , rowind , colind , info )
        if(info/=0) STOP 'info/=0 in unpack_csr_to_coo'
   end subroutine
end module
#endif
