#ifdef CPP_MKL
module m_energy_solve_sparse
use MKL_SPBLAS
use m_tb_types
USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_DOUBLE_COMPLEX,C_PTR,C_F_POINTER
implicit none
private
public Hr_eigvec_sparse_feast, Hr_eigval_sparse_feast
contains
    subroutine Hr_eigvec_sparse_feast(h_par,Hr,eigvec,eigval)
        type(parameters_TB_Hsolve),intent(in)   ::  h_par
        type(SPARSE_MATRIX_T),intent(in)   :: Hr
        complex(8),intent(out),allocatable :: eigvec(:,:)
        real(8),intent(out),allocatable    :: eigval(:)

        !export csr parameters
        integer(C_INT)                     :: nrows,ncols,indexing
        type(C_PTR)                        :: rows_start,rows_end,col_indx
        type(C_PTR)                        :: values
        integer                            :: nnz
        integer,pointer                    :: tmp(:)
        complex(8),pointer                 :: a(:)
        integer,pointer                    :: ia(:),ja(:)

        integer(C_int)              :: stat

        integer                     :: fpm(128)
        real(8),allocatable         :: e(:)
        complex(8),allocatable      :: x(:,:)
        real(8)                     :: emin,emax
        real(8)                     :: epsout
        integer                     :: loop
        integer                     :: m0,m
        real(8),allocatable         :: res(:)
        integer                     :: info


        stat=mkl_sparse_z_export_csr ( Hr , indexing , nrows , ncols , rows_start , rows_end , col_indx , values ) 
        CAll C_F_POINTER(rows_end, tmp,[nrows] )
        nnz=tmp(nrows)-indexing
        CAll C_F_POINTER(values, a,[nnz]) 
        CAll C_F_POINTER(col_indx, ja,[nnz]) 
        CAll C_F_POINTER(rows_start, ia,[h_par%dimH+1])

        Call feastinit(fpm) 
        fpm(1)=1
        fpm(2)=8
        emin=h_par%extE(1)
        emax=h_par%extE(2)
        m0=h_par%estNe
        if(m0==0.or.m0>h_par%dimH) m0=h_par%dimH
        allocate(e(m0),source=0.0d0)
        allocate(x(h_par%dimh,m0),source=cmplx(0.0d0,0.0d0,8))
        allocate(res(h_par%dimH),source=0.0d0)
        call zfeast_hcsrev ( 'F' , h_par%dimH , a , ia , ja , fpm , epsout , loop , emin , emax , m0 , e , x , m , res , info )
        if(info/=0) STOP "sparse diagonalization failed"
        allocate(eigval(m),source=0.0d0)
        eigval=e(1:m)
        allocate(eigvec(h_par%dimH,m),source=cmplx(0.0d0,0.0d0,8))
        eigvec=x(1:h_par%dimH,1:m)

    end subroutine 

    subroutine Hr_eigval_sparse_feast(h_par,Hr,eigval)
        type(parameters_TB_Hsolve),intent(in)   ::  h_par
        type(SPARSE_MATRIX_T),intent(in)   :: Hr
        real(8),intent(out),allocatable    :: eigval(:)

        complex(8),allocatable             :: eigvec(:,:)

        Call Hr_eigvec_sparse_feast(h_par,Hr,eigvec,eigval)


    end subroutine 




end module
#endif
