#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cufft.h>
#include <iostream>
#include <cassert>

__global__ void k_time_m(
    cufftDoubleComplex* h, 
    cufftDoubleComplex* m,
    cufftDoubleComplex* k, 
    int dim_mode,
    int N){

    int i = blockDim.x * blockIdx.x + threadIdx.x;
    int n = (i/dim_mode)*dim_mode;
    if (i < N*dim_mode){
        h[i]=make_cuDoubleComplex(0.0,0.0);
        for (int j=0;j<dim_mode;j++){
            h[i] = cuCadd(h[i] ,cuCmul(k[i*dim_mode+j], m[n+j]));
        }
    }
}    

__global__ void add_cmplx(
    cufftDoubleComplex* arr_sum, 
    cufftDoubleComplex* arr_add,
    int N){

    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < N){
        arr_sum[i]=cuCadd(arr_sum[i] ,arr_add[i]);
    }
}    



extern "C"{
void cuda_fft_calc_h(
    const int& dim_mode,
    int N_rep[3],
    cufftDoubleReal*    &m_n,
    cufftDoubleComplex* &m_f,
    cufftDoubleComplex* &k_f,
    cufftDoubleReal*    &h_n,
    cufftDoubleComplex* &h_f,
    cufftHandle*        &plan_fwd,
    cufftHandle*        &plan_bwd){

    // get magnetic operator in fourier space
    cufftResult result=cufftExecD2Z(*plan_fwd,m_n,m_f);
    assert(result == CUFFT_SUCCESS);

    cudaError_t err=cudaDeviceSynchronize();
    assert(err == cudaSuccess);

    int NK_tot=N_rep[0]*N_rep[1]*N_rep[2];
    int N=NK_tot*dim_mode;

    int threadsPerBlock = std::min(256,N);
    int numBlocks= (N+ threadsPerBlock - 1) / threadsPerBlock;

    k_time_m<<<numBlocks, threadsPerBlock>>>(h_f, m_f, k_f, dim_mode, NK_tot);

    err=cudaDeviceSynchronize();
    assert(err == cudaSuccess);
    result=cufftExecZ2D(*plan_bwd,h_f,h_n);
    assert(result == CUFFT_SUCCESS);
    err=cudaDeviceSynchronize();
    assert(err == cudaSuccess);
}

void cuda_fft_add_cmplx(
    const int N,
    cufftDoubleComplex* &arr_sum,
    cufftDoubleComplex* &arr_add){

    int threadsPerBlock = std::min(256,N);
    int numBlocks= (N+ threadsPerBlock - 1) / threadsPerBlock;

    add_cmplx<<<numBlocks, threadsPerBlock>>>(arr_sum, arr_add,N);

    cudaError_t err=cudaDeviceSynchronize();
    assert(err == cudaSuccess);
}

}
