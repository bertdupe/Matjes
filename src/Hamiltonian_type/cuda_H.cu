#include <iostream>
#include <cuda_runtime.h>
#include <cusparse.h>
#include <cassert>

#include <thrust/device_ptr.h>
#include <thrust/fill.h>

#ifdef CPP_MPI
#include <mpi.h>
#endif
using namespace std;

//__global__ void add_col_mat(
//    double* res, 
//    int* row,
//    int* col,
//    double* val,
//    int col_ind,
//    int ind,
//    double val_vec){
//
//    int i = blockDim.x * blockIdx.x + threadIdx.x;
//    int Nentry=row[col_ind+1]-row[col_ind];
//    int col_ind;
////    printf("Nentry: %d \n", Nentry);
//    if (i < Nentry){
//        int col_ind=row[col_ind]+i;
//        res[col_ind]=res[col_ind] + val[col_ind] * val_vec;
//    }
//}    



extern "C"{

void cuda_create_handle(
    cusparseHandle_t* &handle
    ){
    handle = new cusparseHandle_t;
    cusparseStatus_t stat = cusparseCreate(handle);
    assert (CUSPARSE_STATUS_SUCCESS == stat);
}

void cuda_H_init(
    // initializes cuda csr 
    const int nnz,
    const int Hdim[2],
    int rows[],
    int cols[],
    double arr_in[],
    cusparseSpMatDescr*& spMatDescr,
    cusparseHandle_t*& handle
    ){

    //define some initial helping variables
    cusparseStatus_t stat;
    cudaError_t err;
    size_t buffersize;
    int nrows=Hdim[0];
    int ncols=Hdim[1];

    //sorting the coo-matrix input into the format necessary for CUDA
      //prepare permutation array for sorting on GPU
    stat = cusparseXcoosort_bufferSizeExt(*handle, nrows, ncols, nnz, rows, cols, &buffersize);
    assert (CUSPARSE_STATUS_SUCCESS == stat);
    cudaDeviceSynchronize();
    void *pBuffer = NULL;
    cudaMalloc( &pBuffer, (size_t) sizeof(char)* buffersize);
    int *P = NULL;
    cudaMalloc( (void**)&P, (size_t) sizeof(int)*nnz);
    stat=cusparseCreateIdentityPermutation(*handle, nnz, P);
    assert (CUSPARSE_STATUS_SUCCESS == stat);
    cudaDeviceSynchronize();

      //Prepare actuall data-array for coo on GPU
    int* cooCols;
    int* cooRows;
    err =cudaMalloc( &cooRows, (size_t) sizeof(int)*nnz);
    assert( cudaSuccess == err );
    err =cudaMalloc( &cooCols, (size_t) sizeof(int)*nnz);
    assert( cudaSuccess == err );
    err = cudaMemcpy(cooRows, rows, (size_t) sizeof(int)*nnz , cudaMemcpyHostToDevice);
    assert( cudaSuccess == err );
    err = cudaMemcpy(cooCols, cols, (size_t) sizeof(int)*nnz , cudaMemcpyHostToDevice);
    assert( cudaSuccess == err );
    cudaDeviceSynchronize();

      //Finally sort the cooRows and cooCols
    stat = cusparseXcoosortByRow(*handle, nrows, ncols, nnz, cooRows, cooCols, P, pBuffer);
    assert( cudaSuccess == err );
    cudaDeviceSynchronize();

    //Prepare and get csrRow
    int* csrRow;
    err =cudaMalloc( &csrRow, (size_t) sizeof(int)*(nrows+1));
    assert( cudaSuccess == err );
    cudaDeviceSynchronize();
    stat =  cusparseXcoo2csr(*handle, cooRows, nnz, nrows, csrRow, CUSPARSE_INDEX_BASE_ONE);
    assert( cudaSuccess == err );
    cudaDeviceSynchronize();

    //Move values to GPU
    double* values;
    err =cudaMalloc( &values, (size_t) sizeof(double)*nnz);
    assert( cudaSuccess == err );
    err = cudaMemcpy(values, arr_in, (size_t) sizeof(double)*nnz , cudaMemcpyHostToDevice);
    assert( cudaSuccess == err );


    //Create CSR sparse matrix with Handle
    stat= cusparseCreateCsr(&spMatDescr, nrows, ncols, nnz, csrRow, cooCols, values, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ONE, CUDA_R_64F);
    assert( CUSPARSE_STATUS_SUCCESS == stat);
    cudaDeviceSynchronize();

    //Free stuff no longer necessary
    cudaFree(pBuffer);
    cudaFree(P);
    cudaFree(cooRows);
}


void cuda_H_get_deviceptr(
    cusparseSpMatDescr* &mat,
    void* &csrRow,
    void* &csrCol, 
    void* &csrVal){

    //declare H_in host data
    int64_t rows, cols, nnz;
    cusparseIndexType_t csrRowOffsetsType;
    cusparseIndexType_t csrColIndType;
    cusparseIndexBase_t idxBase;
    cudaDataType        valueType;

    //get data from old sparse reference
    cusparseStatus_t stat = cusparseCsrGet(mat, &rows, &cols, &nnz
                                          ,&csrRow,&csrCol, &csrVal
                                          ,&csrRowOffsetsType, &csrColIndType, &idxBase, &valueType);
    assert( CUSPARSE_STATUS_SUCCESS == stat);
}


void cuda_H_add(
    cusparseSpMatDescr* &mat_1,
    cusparseSpMatDescr* &mat_2,
    cusparseSpMatDescr* &mat_s,
    cusparseHandle_t* &handle){

    //we always only want to add the arrays without prefactor
    const double alpha=1.0;
    const double beta=1.0;

    //declare types for all 3 arrays
    int64_t rows_1, rows_2;
    int64_t cols_1, cols_2;
    int64_t nnz_1, nnz_2, nnz_s;
    int    *csrRow_1=NULL, *csrRow_2=NULL, *csrRow_s=NULL;
    int    *csrCol_1=NULL, *csrCol_2=NULL, *csrCol_s=NULL;
    double *csrVal_1=NULL, *csrVal_2=NULL, *csrVal_s=NULL;
    cusparseStatus_t stat;
    cudaError_t err;

    //set size of arrays to be added
    stat= cusparseSpMatGetSize(mat_1, &rows_1, &cols_1, &nnz_1);
    assert( CUSPARSE_STATUS_SUCCESS == stat );
    stat= cusparseSpMatGetSize(mat_2, &rows_2, &cols_2, &nnz_2);
    assert( CUSPARSE_STATUS_SUCCESS == stat );

    //get pointers to Device data arrays of input matrices
    void *csrRow, *csrCol, *csrVal;
    cuda_H_get_deviceptr(mat_1,csrRow,csrCol,csrVal);
    csrRow_1=(int*) csrRow, csrCol_1=(int*) csrCol, csrVal_1=(double*) csrVal;
    cuda_H_get_deviceptr(mat_2,csrRow,csrCol,csrVal);
    csrRow_2=(int*) csrRow, csrCol_2=(int*) csrCol, csrVal_2=(double*) csrVal;

    //Some neccessary setting apparently
    cusparseSetPointerMode(*handle, CUSPARSE_POINTER_MODE_HOST);

    //Allocate sum row array
    cudaMalloc((void**)&csrRow_s, sizeof(int)*(rows_1+1));

    //Set matrix descriptions for all arrays with index 1
    cusparseMatDescr_t descr_1, descr_2, descr_s;
    stat=cusparseCreateMatDescr(&descr_1);
    assert( CUSPARSE_STATUS_SUCCESS == stat );
    cusparseSetMatIndexBase(descr_1, CUSPARSE_INDEX_BASE_ONE);
    stat=cusparseCreateMatDescr(&descr_2);
    assert( CUSPARSE_STATUS_SUCCESS == stat );
    cusparseSetMatIndexBase(descr_2, CUSPARSE_INDEX_BASE_ONE);
    stat=cusparseCreateMatDescr(&descr_s);
    assert( CUSPARSE_STATUS_SUCCESS == stat );
    cusparseSetMatIndexBase(descr_s, CUSPARSE_INDEX_BASE_ONE);

    //Prepare temporary array (buffer) necessary for addition
    size_t bufferSizeInBytes;
    char *buffer = NULL;
    stat=cusparseDcsrgeam2_bufferSizeExt(*handle, rows_1, cols_1,
        &alpha,
        descr_1, nnz_1,
        csrVal_1, csrRow_1, csrCol_1,
        &beta,
        descr_2, nnz_2,
        csrVal_2, csrRow_2, csrCol_2,
        descr_s,
        csrVal_s, csrRow_s, csrCol_s,
        &bufferSizeInBytes
        );
    assert( CUSPARSE_STATUS_SUCCESS == stat );
    err=cudaMalloc((void**)&buffer, sizeof(char)*bufferSizeInBytes);
    assert( cudaSuccess == err );

    //Get nnz_s 
    int nnzC;
    int *nnzTotalDevHostPtr = &nnzC;
    stat=cusparseXcsrgeam2Nnz(*handle, rows_1, cols_1,
        descr_1, nnz_1, csrRow_1, csrCol_1,
        descr_2, nnz_2, csrRow_2, csrCol_2,
        descr_s, csrRow_s, nnzTotalDevHostPtr,
        buffer);
    assert( CUSPARSE_STATUS_SUCCESS == stat );
    if (NULL != nnzTotalDevHostPtr){
        nnz_s = *nnzTotalDevHostPtr;
    }else{
        int baseC;
        cudaMemcpy(&nnz_s, csrRow_s+rows_1, sizeof(int), cudaMemcpyDeviceToHost);
        cudaMemcpy(&baseC, csrRow_s,        sizeof(int), cudaMemcpyDeviceToHost);
        nnz_s -= baseC;
    }

    //allocate nnz_s dependent arrays (Column and value)
    cudaMalloc((void**)&csrCol_s, sizeof(int)*nnz_s);
    cudaMalloc((void**)&csrVal_s, sizeof(double)*nnz_s);

    //Finally get addition
    stat=cusparseDcsrgeam2(*handle, rows_1, cols_1,
        &alpha,
        descr_1, nnz_1,
        csrVal_1, csrRow_1, csrCol_1,
        &beta,
        descr_2, nnz_2,
        csrVal_2, csrRow_2, csrCol_2,
        descr_s,
        csrVal_s, csrRow_s, csrCol_s,
        buffer);
    assert( CUSPARSE_STATUS_SUCCESS == stat );

    //Create the Matrix handle
    stat= cusparseCreateCsr(&mat_s, rows_1, cols_1, nnz_s, csrRow_s, csrCol_s, csrVal_s, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ONE, CUDA_R_64F);
    assert( CUSPARSE_STATUS_SUCCESS == stat);

    //free temporary allocated variables
    err=cudaFree(buffer);
    assert( cudaSuccess == err );
}

void cuda_H_copy(
    cusparseSpMatDescr* &H_in,
    cusparseSpMatDescr* &H_copy){


    //declare H_in host data
    int64_t rows, cols, nnz;
    cusparseIndexType_t csrRowOffsetsType;
    cusparseIndexType_t csrColIndType;
    cusparseIndexBase_t idxBase;
    cudaDataType        valueType;
    //declare H_in device data
    void *csrRow, *csrCol, *csrVal;


    //get data from old sparse reference
    cusparseStatus_t stat = cusparseCsrGet(H_in, &rows, &cols, &nnz
                                          ,&csrRow,&csrCol, &csrVal
                                          ,&csrRowOffsetsType, &csrColIndType, &idxBase, &valueType);
    assert( CUSPARSE_STATUS_SUCCESS == stat);

    //allocate new data arrays and copy data from previous type
    void *row_new, *col_new, *val_new;
    cudaError_t err;
    err =cudaMalloc( &col_new, (size_t) sizeof(int)*nnz);
    assert( cudaSuccess == err );
    err = cudaMemcpy(col_new, csrCol, (size_t) sizeof(int)*nnz , cudaMemcpyDeviceToDevice);
    assert( cudaSuccess == err );

    err =cudaMalloc( &row_new, (size_t) sizeof(int)*(rows+1));
    assert( cudaSuccess == err );
    err = cudaMemcpy(row_new, csrRow, (size_t) sizeof(int)*(rows+1) , cudaMemcpyDeviceToDevice);
    assert( cudaSuccess == err );

    err =cudaMalloc( &val_new, (size_t) sizeof(double)*nnz);
    assert( cudaSuccess == err );
    err = cudaMemcpy(val_new, csrVal, (size_t) sizeof(double)*nnz , cudaMemcpyDeviceToDevice);
    assert( cudaSuccess == err );

    //create new sparse matrix reference
    stat = cusparseCreateCsr(&H_copy, rows, cols, nnz, row_new, col_new, val_new, csrRowOffsetsType, csrColIndType, idxBase, valueType);
    assert( CUSPARSE_STATUS_SUCCESS == stat);
}


void cuda_H_destroy(
    cusparseSpMatDescr* &mat ){

    //Host data (should get destroyed with cusparseDestroySpMat)
    int64_t rows, cols, nnz;
    cusparseIndexType_t csrRowOffsetsType;
    cusparseIndexType_t csrColIndType;
    cusparseIndexBase_t idxBase;
    cudaDataType        valueType;
    //Device data which must be freed first
    void *csrRow, *csrCol, *csrVal;

    //get data from old sparse reference
    cusparseStatus_t stat = cusparseCsrGet(mat, &rows, &cols, &nnz
                                          ,&csrRow,&csrCol, &csrVal
                                          ,&csrRowOffsetsType, &csrColIndType, &idxBase, &valueType);
    assert( CUSPARSE_STATUS_SUCCESS == stat);

    cudaError_t err;
    err=cudaFree(csrRow);
    assert(err == cudaSuccess);
    err=cudaFree(csrCol);
    assert(err == cudaSuccess);
    err=cudaFree(csrVal);
    assert(err == cudaSuccess);

    stat=cusparseDestroySpMat(mat);

    assert( CUSPARSE_STATUS_SUCCESS == stat );
}


void cuda_H_mult_r_beta(
    cusparseSpMatDescr* &mat,
    cusparseDnVecDescr* &in_vec,
    cusparseDnVecDescr* &out_vec,
    const double        &beta,
    char*               &buffer,
    cusparseHandle_t*   &handle){

    const double alpha= 1.0;    //save in device memory
    //do the actual multiplication
    auto stat = cusparseSpMV(*handle,
                        CUSPARSE_OPERATION_NON_TRANSPOSE,
                        &alpha,
                        mat,
                        in_vec,
                        &beta,
                        out_vec,
                        CUDA_R_64F,
                        CUSPARSE_SPMV_CSR_ALG1,
                        buffer);
#ifdef CPP_DEBUG
    assert (CUSPARSE_STATUS_SUCCESS == stat);
#endif
}

void cuda_H_mult_r(
    cusparseSpMatDescr* &mat,
    cusparseDnVecDescr* &in_vec,
    cusparseDnVecDescr* &out_vec,
    char*               &buffer,
    cusparseHandle_t*   &handle){

    const double alpha= 1.0;    //save in device memory
    const double beta = 0.0;
    //do the actual multiplication
    auto stat = cusparseSpMV(*handle,
                        CUSPARSE_OPERATION_NON_TRANSPOSE,
                        &alpha,
                        mat,
                        in_vec,
                        &beta,
                        out_vec,
                        CUDA_R_64F,
                        CUSPARSE_SPMV_CSR_ALG1,
                        buffer);
#ifdef CPP_DEBUG
    assert (CUSPARSE_STATUS_SUCCESS == stat);
#endif
}


void cuda_H_mult_l(
    cusparseSpMatDescr* &mat,
    cusparseDnVecDescr* &in_vec,
    cusparseDnVecDescr* &out_vec,
    char*               &buffer,
    cusparseHandle_t*   &handle){

    const double alpha= 1.0;    //save in device memory
    const double beta = 0.0;
    //do the actual multiplication
    auto stat = cusparseSpMV(*handle,
                        CUSPARSE_OPERATION_TRANSPOSE,
                        &alpha,
                        mat,
                        in_vec,
                        &beta,
                        out_vec,
                        CUDA_R_64F,
                        CUSPARSE_SPMV_CSR_ALG1,
                        buffer);
#ifdef CPP_DEBUG
    assert (CUSPARSE_STATUS_SUCCESS == stat);
#endif
}

void cuda_H_mult_l_beta(
    cusparseSpMatDescr* &mat,
    cusparseDnVecDescr* &in_vec,
    cusparseDnVecDescr* &out_vec,
    const double        &beta,
    char*               &buffer,
    cusparseHandle_t*   &handle){

    const double alpha= 1.0;    //save in device memory
//    const double beta = 0.0;
    //do the actual multiplication
    auto stat = cusparseSpMV(*handle,
                        CUSPARSE_OPERATION_TRANSPOSE,
                        &alpha,
                        mat,
                        in_vec,
                        &beta,
                        out_vec,
                        CUDA_R_64F,
                        CUSPARSE_SPMV_CSR_ALG1,
                        buffer);
#ifdef CPP_DEBUG
    assert (CUSPARSE_STATUS_SUCCESS == stat);
#endif
}


void cuda_free_buffer(
    char* &buffer
    ){

    auto err=cudaFree(buffer);
    assert( cudaSuccess == err);
}


void cuda_set_buffer(
    char* &buffer,
    cusparseSpMatDescr* &mat,
    bool transpose,
    cusparseDnVecDescr* &in_vec,
    cusparseDnVecDescr* &out_vec,
    cusparseHandle_t* &handle){

    cusparseStatus_t stat;
    cudaError_t err;

    //get device buffer
    size_t bufferSize=0;
    const double alpha= 1.0;
    const double beta = 0.0;
    cusparseOperation_t operation;
    if(transpose){
        operation=CUSPARSE_OPERATION_TRANSPOSE;
    }else{
        operation=CUSPARSE_OPERATION_NON_TRANSPOSE;
    }

    stat =cusparseSpMV_bufferSize(*handle,
                                  operation,
                                  &alpha,
                                  mat,
                                  in_vec,
                                  &beta,
                                  out_vec,
                                  CUDA_R_64F,
                                  CUSPARSE_SPMV_CSR_ALG1,
                                  &bufferSize);
    assert (CUSPARSE_STATUS_SUCCESS == stat);
    err=cudaMalloc( (void**) &buffer, (size_t) sizeof(char)* bufferSize);
    assert( cudaSuccess == err);
}

//void cuda_H_mult_mat_disc_disc(
//    //THIS ROUTINE IS NOT FINISHED AND PROBABLY NEVER WILL BE SINCE DIRECT EVALUATION OF THE SPARSE MATRICES IS FASTER DIRECTLY ON THE CPU
//    cusparseSpMatDescr* &mat,
//    cusparseDnVecDescr* &in_vec,
//    cusparseDnVecDescr* &out_vec,
//    int N_in,
//    int ind_in[],
//    double vec_in[],
//    int& N_out,  // in: size of ind_out , out: size of relevant indices in ind_out
//    int ind_out[],
//    double vec_out[]){
//
//    cudaError_t err;
//
//    //host data
//    int64_t rows, cols, nnz;
//    cusparseIndexType_t csrRowOffsetsType;
//    cusparseIndexType_t csrColIndType;
//    cusparseIndexBase_t idxBase;
//    cudaDataType        valueType;
//    //Device data
//    void *csrRow_v, *csrCol_v, *csrVal_v;
//    int *csrRow, *csrCol;
//    double *csrVal;
//
//    //get data from old sparse reference
//    cusparseStatus_t stat = cusparseCsrGet(mat, &rows, &cols, &nnz
//                                          ,&csrRow_v,&csrCol_v, &csrVal_v
//                                          ,&csrRowOffsetsType, &csrColIndType, &idxBase, &valueType);
//    assert( CUSPARSE_STATUS_SUCCESS == stat );
//
//    csrRow=(int*) csrRow_v, csrCol=(int*) csrCol_v, csrVal=(double*) csrVal_v;
//
//    //get data from input vector (only used for storage)
//    double* vec_val_dev;
//    stat=cusparseDnVecGetValues(in_vec, (void**) &vec_val_dev);
//    assert( CUSPARSE_STATUS_SUCCESS == stat );
//    thrust::device_ptr<double> dev_ptr(vec_val_dev);
//    thrust::fill(dev_ptr, dev_ptr + cols, 0.0);
//
////    //set in vector
////    double *in_val;
////    int *in_ind;
////    err =cudaMalloc( &in_val,(size_t) sizeof(double)*N_in);
////    assert( cudaSuccess == err );
////    err = cudaMemcpy(in_val, vec_in, (size_t) sizeof(int)*N_in , cudaMemcpyHostToDevice);
////    assert( cudaSuccess == err );
////
////    err =cudaMalloc( &in_ind,(size_t) sizeof(int)*N_in);
////    assert( cudaSuccess == err );
////    err = cudaMemcpy(in_ind, ind_in, (size_t) sizeof(int)*N_in , cudaMemcpyHostToDevice);
////    assert( cudaSuccess == err );
//
//
//    int N=cols;
//    int threadsPerBlock = std::min(256,N);
//    int numBlocks= (N+ threadsPerBlock - 1) / threadsPerBlock;
//
////    for (int i =0; i<N_in; ++i){
////        add_col_mat<<<numBlocks, threadsPerBlock>>>(vec_val_dev,csrRow, csrCol, csrVal, i, ind_in[i], vec_in[i]);
////    }
//}

}
