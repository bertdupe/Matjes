#ifdef CPP_MKL_SPBLAS
module m_energy_solve_sparse
use MKL_SPBLAS
use mkl_spblas_util, only: unpack_csr 
use m_tb_types
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
        complex(8),pointer                 :: a(:)
        integer,pointer                    :: ia(:),ja(:)
        integer                            :: nnz

        integer                     :: fpm(128)
        real(8),allocatable         :: e(:)
        complex(8),allocatable      :: x(:,:)
        real(8)                     :: emin,emax
        real(8)                     :: epsout
        integer                     :: loop
        integer                     :: m0,m
        real(8),allocatable         :: res(:)
        integer                     :: info


        Call unpack_csr(h_par%dimH,Hr,nnz,ia,ja,a)
        Call feastinit(fpm) 
        fpm(1)=1
        fpm(2)=8
        fpm(3)=-nint(log10(h_par%diag_acc))
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
