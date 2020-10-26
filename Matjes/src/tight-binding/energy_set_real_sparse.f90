#ifdef CPP_MKL

module m_energy_set_real_sparse
use m_energy_commons, only : energy
use m_basic_types, only : vec_point
use m_tb_types
use MKL_SPBLAS
use mkl_spblas_util, only: unpack_csr_to_coo ,unpack_csr
USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_DOUBLE_COMPLEX

private
public set_Hr_sparse_nc!set_Jsd!,Hr_eigval,Hr_eigvec


contains

    subroutine set_Hr_sparse_nc(h_par,mode_mag,H_r)
        use m_tb_params, only : TB_params
        type(parameters_TB_Hsolve),intent(in)   ::  h_par
        type(SPARSE_MATRIX_T),intent(out)       ::  H_r
        type(vec_point),intent(in)              ::  mode_mag(:)

        type(SPARSE_MATRIX_T)        :: H_ee
        type(SPARSE_MATRIX_T)        :: H_jsd
        integer(C_int)               :: stat
        type(MATRIX_DESCR)           :: desc



        Call set_Hr_ee(h_par,H_ee)
        if(any(TB_params%io_H%Jsd /= 0.0d0))then
            Call set_Hr_Jsd(h_par,mode_mag,TB_params%io_H%Jsd,H_jsd)
            stat=MKL_SPARSE_Z_ADD(SPARSE_OPERATION_NON_TRANSPOSE,H_ee,cmplx(1.0,0.0,C_DOUBLE_COMPLEX),H_jsd,H_r)
            if(stat /= 0) STOP "error adding sparse H_ee and H_jsd"
            stat=MKL_SPARSE_DESTROY(H_ee)
            stat=MKL_SPARSE_DESTROY(H_jsd)
        else
            H_r=H_ee
        endif
        stat = mkl_sparse_order ( H_r )
        !not sure if that enhances things
        desc%Type=SPARSE_MATRIX_TYPE_HERMITIAN
        desc%mode=SPARSE_FILL_MODE_UPPER
        desc%diag=SPARSE_DIAG_NON_UNIT
        !commented out because I'm not sure if it helps and the interface changes between mkl versions
        !stat = mkl_sparse_set_lu_smoother_hint ( H_r , SPARSE_OPERATION_NON_TRANSPOSE , desc , 100000 )
        stat = mkl_sparse_optimize ( H_r )

   end subroutine 

   subroutine set_Hr_ee(h_par,H_csr)
        type(parameters_TB_Hsolve),intent(in)     ::  h_par
        type(SPARSE_MATRIX_T),intent(out)         ::  H_csr
        !external SPARSE_MATRIX_T
        integer                 :: dimH
        integer                 :: TB_ext(2)
        integer                 :: N_persite

        integer                 :: N_neighbours,N_cells,dim_mode
        real(8)                 :: E
        integer                 :: i,j,k,ii,i2,k2

!sparse Hamiltonian
        integer(C_INT)                          :: nnz
        complex(C_DOUBLE_COMPLEX),allocatable   :: val(:)
        integer(C_INT),allocatable              :: rowind(:),colind(:)
        integer(C_INT)                          :: cols,rows


        type(SPARSE_MATRIX_T)   ::  H_coo
        integer(C_INT)          ::  stat

        dimH=h_par%dimH
        TB_ext=h_par%pos_ext

        N_neighbours = size( energy%line, 1 )
        N_cells = size( energy%line, 2 )
        dim_mode=Tb_ext(2)-Tb_ext(1)+1

        N_persite=0
        do i=1,N_neighbours
            k = energy%line(i, 1) 
            N_persite=N_persite+count(energy%value(i,k)%order_op(1)%Op_loc(TB_ext(1):TB_ext(2),TB_ext(1):TB_ext(2)) /= 0.0d0)
        enddo
        nnz=N_persite*N_cells

        !set the hamiltonian in coordinate matrix storage
        allocate(rowind(nnz),colind(nnz),source=0)
        allocate(val(nnz),source=cmplx(0.0d0,0.0d0,8))
        ii=1
        do i=1, N_cells
            do j=1, N_neighbours
                k = energy%line(j, i)
                do i2=TB_ext(1),TB_ext(2)
                    do k2=TB_ext(1),TB_ext(2)
                        E=energy%value(j,k)%order_op(1)%Op_loc(i2,k2)
                        if(E/=0.0d0)then
                            rowind(ii)=dim_mode*(i-1)+i2-TB_ext(1)+1
                            colind(ii)=dim_mode*(k-1)+k2-TB_ext(1)+1
                            val(ii)=cmplx(E,0.0d0,8)
                            ii=ii+1
                        endif
                    enddo
                enddo
            enddo
        enddo
        cols=dimH;rows=dimH
        stat = mkl_sparse_z_create_coo ( H_coo , SPARSE_INDEX_BASE_ONE , rows , cols , nnz , rowind , colind , val)
        if(stat /= 0) STOP "error creating H_coo in sparse H_ee"
        stat = MKL_SPARSE_CONVERT_CSR(H_coo,SPARSE_OPERATION_NON_TRANSPOSE,H_csr)
        if(stat /= 0) STOP "error setting H_ee sparse to CSR"
        stat=MKL_SPARSE_DESTROY(H_coo)
    end subroutine

    subroutine set_Hr_Jsd(h_par,mode_mag,Jsd,H_csr)
        !use m_eigen_interface, only : eigen_set_H_e_jsd
        !adds the Jsd coupling to a local real-space Hamiltonian Hr
        type(parameters_TB_Hsolve),intent(in)     ::  h_par
        type(SPARSE_MATRIX_T),intent(out)         ::  H_csr
        real(8),intent(in)           ::  Jsd(:)
        type(vec_point),intent(in)   ::  mode_mag(:)

        integer                 :: nnz
        complex(8),allocatable  :: val(:)
        integer,allocatable     :: rowind(:),colind(:)

        complex(8)              :: jsdmat(2,2)
        integer                 :: matind(2)

        integer                 ::  i,j
        integer                 ::  i1,i2,ii

        type(SPARSE_MATRIX_T)   ::  H_coo
        integer(C_INT)          ::  stat

        nnz=h_par%ncell*h_par%norb*4
        allocate(rowind(nnz),colind(nnz),source=0)
        allocate(val(nnz),source=cmplx(0.0d0,0.0d0,8))

        ii=1
        do i=1,h_par%ncell
            do j=1,h_par%norb
                matind=[2*j-1,2*j]+(i-1)*h_par%norb*h_par%nspin
                jsdmat(1,1)=Jsd(j)*cmplx( mode_mag(i)%w(3), 0.0d0           ,8)
                jsdmat(2,1)=Jsd(j)*cmplx( mode_mag(i)%w(1), mode_mag(i)%w(2),8)
                jsdmat(1,2)=Jsd(j)*cmplx( mode_mag(i)%w(1),-mode_mag(i)%w(2),8)
                jsdmat(2,2)=Jsd(j)*cmplx(-mode_mag(i)%w(3), 0.0d0           ,8)
                do i1=1,2
                    do i2=1,2
                        rowind(ii)=matind(i1)
                        colind(ii)=matind(i2)
                        val(ii)=jsdmat(i1,i2)
                        ii=ii+1
                    enddo
                enddo
            enddo
        enddo

        stat = mkl_sparse_z_create_coo ( H_coo , SPARSE_INDEX_BASE_ONE , h_par%dimH , h_par%dimH , nnz , rowind , colind , val)
        if(stat /= 0) STOP "error creating H_coo for sparse H_Jsd"
        stat = MKL_SPARSE_CONVERT_CSR(H_coo,SPARSE_OPERATION_NON_TRANSPOSE,H_csr)
        if(stat /= 0) STOP "error setting H_Jsd sparse to CSR"
        stat=MKL_SPARSE_DESTROY(H_coo)
    end subroutine 
end module
#endif
