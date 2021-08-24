include 'mkl_spblas.f90'

program test_mkl_spblas
    use MKL_SPBLAS, only: mkl_sparse_d_create_coo,SPARSE_MATRIX_T,SPARSE_STATUS_SUCCESS,SPARSE_INDEX_BASE_ONE,mkl_sparse_destroy
    USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_DOUBLE,C_INT

    type(SPARSE_MATRIX_T)       :: H
    integer(C_int)              :: stat
    real(C_DOUBLE)              :: val(3)
    integer                     :: rowind(3),colind(3)

    rowind=[1,1,2]
    colind=[1,2,2]
    val=[1.0d0,2.0d0,3.0d0]

    stat=mkl_sparse_d_create_coo(H, SPARSE_INDEX_BASE_ONE , 2 , 2 , 3 , rowind , colind , val)
    if(stat/=SPARSE_STATUS_SUCCESS) STOP "failed to initialize MKL_SPBLAS matrix"
    stat=mkl_sparse_destroy(H)
    if(stat/=SPARSE_STATUS_SUCCESS) STOP "failed to destory MKL_SPBLAS matrix"
end program    

