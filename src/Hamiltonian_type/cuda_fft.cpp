#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cufft.h>
#include <iostream>
#include <cassert>

extern "C"{

void cuda_fft_init(
    const int& dim_mode,
    const int N_rep[3],
    cufftDoubleReal*    const& m_n,
    cufftDoubleComplex* const& m_f,
    cufftDoubleReal*    const& h_n,
    cufftDoubleComplex* const& h_f){

    int NK_tot=N_rep[0]*N_rep[1]*N_rep[2];
 
    cudaError_t err;
    err=cudaMalloc((void**) &m_n, sizeof(cufftDoubleReal)   *NK_tot*dim_mode);
    assert(err == cudaSuccess);
    
    err=cudaMalloc((void**) &m_f, sizeof(cufftDoubleComplex)*NK_tot*dim_mode);
    assert(err == cudaSuccess);

    err=cudaMalloc((void**) &h_n, sizeof(cufftDoubleReal)   *NK_tot*dim_mode);
    assert(err == cudaSuccess);

    err=cudaMalloc((void**) &h_f, sizeof(cufftDoubleComplex)*NK_tot*dim_mode);
    assert(err == cudaSuccess);
}


void cuda_fft_get_real(
    const int size,
    cufftDoubleReal* const& arr_device,
    double* arr_out){

    cudaError_t err;
#ifdef CPP_DEBUG
    cudaPointerAttributes attributes;
    err=cudaPointerGetAttributes (&attributes, arr_device);
    assert(err == cudaSuccess);
    assert(attributes.type == cudaMemoryTypeDevice);
#endif
    err= cudaMemcpy(arr_out, arr_device,(size_t) size * sizeof(cufftDoubleReal),cudaMemcpyDeviceToHost);
    assert(err == cudaSuccess);
}


void cuda_fft_set_real(
    const int& size,
    const double* arr_in,
    cufftDoubleReal* const& arr_device){

#ifdef CPP_DEBUG
    cudaPointerAttributes attributes;
    cudaError_t err;
    err=cudaPointerGetAttributes (&attributes, arr_device);
    assert(err == cudaSuccess);
    assert(attributes.type == cudaMemoryTypeDevice);
#endif

    err= cudaMemcpy(arr_device, arr_in,(size_t) size * sizeof(cufftDoubleReal),cudaMemcpyHostToDevice);
    assert(err == cudaSuccess);
}


void set_cuda_fft_plan(
    const int& dim_mode,
    int N_rep[3],
    cufftHandle*& plan_fwd,
    cufftHandle*& plan_bwd){
    plan_fwd=new cufftHandle;
    plan_bwd=new cufftHandle;
    cufftResult result;

    int N_rep_r[3] {N_rep[2],N_rep[1],N_rep[0]};    //fortran vs c reverses order
    int embed_r[3] {N_rep_r[0],N_rep_r[1],N_rep_r[2]};
    int embed_c[3] {N_rep_r[0],N_rep_r[1],(N_rep_r[2]/2)+1};

    result=cufftPlanMany(plan_fwd, 3, N_rep_r, 
                         embed_r, dim_mode, 1,
                         embed_c, dim_mode, 1, 
                         CUFFT_D2Z , dim_mode);
    assert(result == CUFFT_SUCCESS);

    result=cufftPlanMany(plan_bwd, 3, N_rep_r, 
                         embed_c, dim_mode, 1,
                         embed_r, dim_mode, 1, 
                         CUFFT_Z2D , dim_mode);
    assert(result == CUFFT_SUCCESS);
}

void cuda_fft_set_operator(
    const int& dim_mode,
    int N_rep[3],
    //const double* arr_n,
    const double* arr_n,
    cufftDoubleComplex* const& op_f){

    //helpfull values
    int NK_tot=N_rep[0]*N_rep[1]*N_rep[2];
    int N=NK_tot*dim_mode*dim_mode;

    //get operator in normal space on device
    cufftDoubleReal* op_n;
    cudaError_t err;
    err=cudaMalloc((void**) &op_n, (size_t) sizeof(cufftDoubleReal) * N);
    assert(err == cudaSuccess);
    err= cudaMemcpy(op_n, arr_n, (size_t) sizeof(cufftDoubleReal) * N,cudaMemcpyHostToDevice);
    assert(err == cudaSuccess);

    //get transformation normal to fourier space
    cufftHandle plan;
    int N_rep_r[3] {N_rep[2],N_rep[1],N_rep[0]};
    int embed_r[3] {N_rep_r[0],N_rep_r[1],N_rep_r[2]};
    int embed_c[3] {N_rep_r[0],N_rep_r[1],(N_rep_r[2]/2)+1};

    cufftResult result=cufftPlanMany(&plan, 3,N_rep_r, 
                        embed_r, dim_mode*dim_mode, 1,
                        embed_c, dim_mode*dim_mode, 1,
                        CUFFT_D2Z , dim_mode*dim_mode);

    assert(result == CUFFT_SUCCESS);

    //get fourier transformed output operator on device
    err=cudaMalloc((void**) &op_f, (size_t) sizeof(cufftDoubleComplex) * N);
    assert(err == cudaSuccess);

    //do the actual transformation
    result=cufftExecD2Z(plan,op_n,op_f);
    assert(result == CUFFT_SUCCESS);
    err=cudaDeviceSynchronize();
    assert(err == cudaSuccess);
    
    //destroy temporary arrays
    result=cufftDestroy(plan);
    assert(result == CUFFT_SUCCESS);
    err=cudaFree(op_n);
    assert(err == cudaSuccess);
}

void cuda_fft_copy_cmplx(
    const int N,
    cufftDoubleComplex* &arr_in,
    cufftDoubleComplex* &arr_out){

    cudaError_t err;
    err=cudaMalloc((void**) &arr_out, (size_t) sizeof(cufftDoubleComplex) * N);
    assert(err == cudaSuccess);
    err= cudaMemcpy(arr_out, arr_in,(size_t) N * sizeof(cufftDoubleComplex),cudaMemcpyDeviceToDevice);
    assert(err == cudaSuccess);
}

void cuda_fft_destroy_plan(
    cufftHandle*& plan){

    cufftResult result=cufftDestroy(*plan);
    assert(result == CUFFT_SUCCESS);
    delete(plan);
    plan=NULL;
}

void cuda_fft_destroy_real(
    cufftDoubleReal* & arr){

    cudaError_t err= cudaFree(arr);
    assert(err == cudaSuccess);
    arr=NULL;
}

void cuda_fft_destroy_cmplx(
    cufftDoubleComplex* & arr){

    cudaError_t err= cudaFree(arr);
    assert(err == cudaSuccess);
    arr=NULL;
}


}

