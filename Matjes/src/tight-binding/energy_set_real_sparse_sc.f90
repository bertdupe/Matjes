#ifdef CPP_MKL_SPBLAS

module m_energy_set_real_sparse_sc
use m_energy_commons, only : energy
use m_tb_types
use MKL_SPBLAS
use m_derived_types, only: lattice
use mkl_spblas_util, only: unpack_csr 
use  m_energy_set_real_sparse, only: set_Hr_sparse_nc
USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_DOUBLE_COMPLEX,C_PTR,C_F_POINTER

private
public set_Hr_sparse_sc

contains

    subroutine set_Hr_sparse_sc(lat,h_par,h_io,mode_mag,H_r)
        use m_tb_params, only : TB_params
        type(lattice),intent(in)                :: lat
        type(parameters_TB_Hsolve),intent(in)   ::  h_par
        type(parameters_TB_IO_H),intent(in)     :: h_io
        type(SPARSE_MATRIX_T),intent(out)       ::  H_r
        real(8),intent(in)                       ::  mode_mag(:,:)

        type(parameters_TB_Hsolve)              ::  h_par_nc
        type(SPARSE_MATRIX_T)        :: H_nc
        type(SPARSE_MATRIX_T)        :: H_delta
        type(SPARSE_MATRIX_T)        :: H_double

        integer(C_int)               :: stat
        !type(MATRIX_DESCR)           :: desc


        !get Hamiltonian without superconductivity
        h_par_nc=h_par
        h_par_nc%nsc=1
        Call h_par_nc%upd()
        Call set_Hr_sparse_nc(lat,h_par_nc,h_io, mode_mag,H_nc)
        Call get_Hr_double_sc(h_par,h_par_nc,H_nc,H_double)
        Call get_delta(TB_params%io_H%delta,h_par,H_delta)

        stat=MKL_SPARSE_Z_ADD(SPARSE_OPERATION_NON_TRANSPOSE,H_double,cmplx(1.0,0.0,C_DOUBLE_COMPLEX),H_delta,H_r)
        if(stat/=0) STOP 'failed to add H_double and H_delta'
        stat=MKL_SPARSE_DESTROY(H_double)
        stat=MKL_SPARSE_DESTROY(H_delta)
        stat = mkl_sparse_order ( H_r )
   end subroutine 


   subroutine get_Hr_double_sc(h_par,h_par_nc,H_nc,H_double)
        !super clumsy way to double the Hamiltonian adding the creators/destructors for SC going memory expensively through the coo representation
        type(parameters_TB_Hsolve),intent(in)   ::  h_par,h_par_nc
        type(SPARSE_MATRIX_T),intent(inout)     ::  H_nc
        type(SPARSE_MATRIX_T),intent(out)       ::  H_double


        !export csr parameters
        integer                            :: nnz
        complex(C_DOUBLE_COMPLEX),pointer  :: acsr(:)
        integer(C_INT),pointer             :: ia(:),ja(:)

        !convert to coordinate representation
        integer                            :: job(6)
        complex(8),allocatable             :: acoo(:)
        integer,allocatable                :: rowind(:),colind(:)
        type(SPARSE_MATRIX_T)              :: H_coo

        integer(C_INT)                     :: stat
        integer                            :: info
        
        !get csr matrix out of SPARSE_MATRIX_T
        Call unpack_csr(h_par_nc%dimH,H_nc,nnz,ia,ja,acsr)

        !create coo sparse matrix and add negative branch in (2,2) quadrant
        job=[0,1,1,0,nnz*2,3]
        allocate(acoo(nnz*2),source=cmplx(9.0d0,0.0d0,8))
        allocate(rowind(nnz*2),source=0)
        allocate(colind(nnz*2),source=0)
        Call mkl_zcsrcoo ( job , h_par_nc%dimH , acsr , ja , ia , nnz , acoo , rowind , colind , info )
        rowind(nnz+1:2*nnz)=rowind(1:nnz)+h_par%dimH/2
        colind(nnz+1:2*nnz)=colind(1:nnz)+h_par%dimH/2
        acoo(nnz+1:2*nnz)=-acoo(1:nnz)
        nnz=nnz*2

        !create again SPARSE_MATRIX_T type and transform intenally to CSR
        stat = mkl_sparse_z_create_coo ( H_coo , SPARSE_INDEX_BASE_ONE , h_par%dimH , h_par%dimH , nnz , rowind , colind , acoo)
        if(stat /= 0) STOP "error creating H_coo in sparse H_ee"
        stat = MKL_SPARSE_CONVERT_CSR(H_coo,SPARSE_OPERATION_NON_TRANSPOSE,H_double)
        if(stat /= 0) STOP "error setting H_ee sparse to CSR"
        stat=MKL_SPARSE_DESTROY(H_coo)
   end subroutine

   subroutine get_delta(delta,h_par,H_delta)
        complex(8),intent(in)                   ::  delta(:)
        type(parameters_TB_Hsolve),intent(in)   ::  h_par
        type(SPARSE_MATRIX_T),intent(out)       ::  H_delta

        !matrix in coo
        integer                         :: nnz
        complex(8),allocatable          :: acoo(:)
        integer,allocatable             :: rowind(:),colind(:)
        type(SPARSE_MATRIX_T)           :: H_coo

        integer             ::  i_cell,i_orb
        integer             ::  i_up,i_dn,i_up_dg,i_dn_dg
        integer             ::  ii
        integer(C_INT)      :: stat


        nnz=h_par%dimH
        allocate(rowind(nnz),source=0)
        allocate(colind(nnz),source=0)
        allocate(acoo(nnz),source=cmplx(0.0d0,0.0d0,8))
        ii=1
        do i_cell=1,h_par%ncell
            do i_orb=1,h_par%norb
                i_up=2*(i_cell-1)*i_orb+2*i_orb-1
                i_dn=i_up+1
                i_up_dg=i_up+h_par%dimH/2
                i_dn_dg=i_dn+h_par%dimH/2
                !Hr(i_up_dg,i_dn)=Hr(i_up_dg,i_dn)-conjg(delta(i_orb))
                colind(ii)=i_up_dg
                rowind(ii)=i_dn
                acoo(ii)=-conjg(delta(i_orb))
                ii=ii+1
                !hr(i_dn_dg,i_up)=Hr(i_dn_dg,i_up)+conjg(delta(i_orb))
                colind(ii)=i_dn_dg
                rowind(ii)=i_up
                acoo(ii)=+conjg(delta(i_orb))
                ii=ii+1
                !Hr(i_up,i_dn_dg)=Hr(i_up,i_dn_dg)+delta(i_orb)
                colind(ii)=i_up
                rowind(ii)=i_dn_dg
                acoo(ii)=delta(i_orb)
                ii=ii+1
                !Hr(i_dn,i_up_dg)=Hr(i_dn,i_up_dg)-delta(i_orb)
                colind(ii)=i_dn
                rowind(ii)=i_up_dg
                acoo(ii)=-delta(i_orb)
                ii=ii+1
            enddo
        enddo
        stat = mkl_sparse_z_create_coo ( H_coo , SPARSE_INDEX_BASE_ONE , h_par%dimH , h_par%dimH , nnz , rowind , colind , acoo)
        if(stat /= 0) STOP "error creating H_coo in sparse H_ee"
        stat = MKL_SPARSE_CONVERT_CSR(H_coo,SPARSE_OPERATION_NON_TRANSPOSE,H_delta)
        if(stat /= 0) STOP "error setting H_ee sparse to CSR"
        stat=MKL_SPARSE_DESTROY(H_coo)
   end subroutine

end module
#endif
