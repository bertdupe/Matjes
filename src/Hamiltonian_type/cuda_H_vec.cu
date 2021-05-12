#include <iostream>
#include <cuda_runtime.h>
#include <cusparse.h>
#include <cassert>
#ifdef CPP_MPI
#include <mpi.h>
#endif
using namespace std;

extern "C"{


void cuda_fvec_alloccopy(
    cusparseDnVecDescr* &vec_in, 
    cusparseDnVecDescr* &vec_out){

    double *val_in, *val_out;
    //get references
    cudaDataType dat_type;
    int64_t size;
    auto stat=cusparseDnVecGet(vec_in,&size, (void**) &val_in,&dat_type);
    assert (CUSPARSE_STATUS_SUCCESS == stat);

    //allocate and copy new device memory 
    auto err =cudaMalloc( (void**) &val_out, (size_t) sizeof(double)*size);
    assert( cudaSuccess == err );
    err = cudaMemcpy(val_out, val_in, (size_t) sizeof(double)*size , cudaMemcpyDeviceToDevice);
    assert( cudaSuccess == err);
    //initialize type
    stat= cusparseCreateDnVec(&vec_out,  (int64_t) size, val_out,  CUDA_R_64F);
    assert (CUSPARSE_STATUS_SUCCESS == stat);
    }

void cuda_fvec_init(
    cusparseDnVecDescr* &vec, 
    int size){

    double *val;
    //allocate device memory 
    auto err =cudaMalloc( (void**) &val, (size_t) sizeof(double)*size);
    assert( cudaSuccess == err );
    //initialize type
    auto stat= cusparseCreateDnVec(&vec,  (int64_t) size, val,  CUDA_R_64F);
    assert (CUSPARSE_STATUS_SUCCESS == stat);
    }

void cuda_fvec_destroy(
    cusparseDnVecDescr* &vec){

    double *val;
    //free device memory
    auto stat=cusparseDnVecGetValues(vec, (void**) &val);
    assert (CUSPARSE_STATUS_SUCCESS == stat);
    auto err=cudaFree(val);
    assert( cudaSuccess == err);
    
    //free host memory
    stat=cusparseDestroyDnVec(vec);
    assert (CUSPARSE_STATUS_SUCCESS == stat);
}

void cuda_fvec_set(
    cusparseDnVecDescr* &vec,
    const double arr_in[]){

    //get references
    double *val;
    cudaDataType dat_type;
    int64_t size;
    auto stat=cusparseDnVecGet(vec,&size, (void**) &val,&dat_type);
    assert (CUSPARSE_STATUS_SUCCESS == stat);
    //copy array
    auto err = cudaMemcpy(val, arr_in, (size_t) sizeof(double)*size , cudaMemcpyHostToDevice);
    assert( cudaSuccess == err);
}

void cuda_fvec_get(
    cusparseDnVecDescr* &vec,
    double arr_out[]){

    //get references
    double *val;
    cudaDataType dat_type;
    int64_t size;
    auto stat=cusparseDnVecGet(vec,&size, (void**) &val,&dat_type);
    assert (CUSPARSE_STATUS_SUCCESS == stat);
    //copy array
    auto err = cudaMemcpy(arr_out, val, (size_t) sizeof(double)*size , cudaMemcpyDeviceToHost);
    assert( cudaSuccess == err);
}
}
