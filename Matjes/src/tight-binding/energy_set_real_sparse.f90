#ifdef CPP_MKL_SPBLAS

module m_energy_set_real_sparse
use m_tb_types
use MKL_SPBLAS
use m_derived_types, only: lattice
use mkl_spblas_util, only: unpack_csr_to_coo ,unpack_csr
use m_get_table_nn,only :get_table_nn
USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_DOUBLE_COMPLEX

private
public set_Hr_sparse_nc!set_Jsd!,Hr_eigval,Hr_eigvec


contains

    subroutine set_Hr_sparse_nc(lat,h_par,h_io,mode_mag,H_r)
        type(lattice),intent(in)                :: lat
        type(parameters_TB_Hsolve),intent(in)   :: h_par
        type(parameters_TB_IO_H),intent(in)     :: h_io
        type(SPARSE_MATRIX_T),intent(out)       ::  H_r
        real(8),intent(in)                      ::  mode_mag(:,:)
        !local
        type(SPARSE_MATRIX_T)        :: H_ee
        type(SPARSE_MATRIX_T)        :: H_jsd
        integer(C_int)               :: stat
        type(MATRIX_DESCR)           :: desc

        Call set_Hr_ee(lat,h_par,h_io,H_ee)
        if(any(h_io%Jsd /= 0.0d0))then
            Call set_Hr_Jsd(h_par,mode_mag,h_io%Jsd,H_jsd)
            stat=MKL_SPARSE_Z_ADD(SPARSE_OPERATION_NON_TRANSPOSE,H_ee,cmplx(1.0,0.0,C_DOUBLE_COMPLEX),H_jsd,H_r)
            if(stat /= 0) ERROR STOP "error adding sparse H_ee and H_jsd"
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

   subroutine set_Hr_ee(lat,h_par,h_io,H_csr)
        type(lattice),intent(in)                    :: lat
        type(parameters_TB_IO_H),intent(in)     :: h_io
        type(parameters_TB_Hsolve),intent(in)     ::  h_par
        type(SPARSE_MATRIX_T),intent(out)         ::  H_csr
        !external SPARSE_MATRIX_T
        integer                 :: dimH
        integer                 :: N_persite

        integer                 :: N_neighbours,Ncell,dim_mode
        integer                 :: i,ii
        integer, allocatable    :: indexNN(:),tableNN(:,:,:,:,:,:)
        integer                 :: i_cell,i_orb,i_s,i_neigh,i_vois
        integer                 :: ind_orb
        integer                 :: offset
        integer                 :: ilat_1(3),ilat_2(3),i_other
        complex(8)              :: H_entry

!sparse Hamiltonian
        integer(C_INT)                          :: nnz
        complex(C_DOUBLE_COMPLEX),allocatable   :: val(:)
        integer(C_INT),allocatable              :: rowind(:),colind(:)
        integer(C_INT)                          :: cols,rows


        type(SPARSE_MATRIX_T)   ::  H_coo
        integer(C_INT)          ::  stat


        N_neighbours = size(h_io%hopping,3)
        Call get_table_nn(lat,N_neighbours,indexNN,tableNN)

        dimH=h_par%dimH
        Ncell = lat%Ncell
        dim_mode=h_par%norb*h_par%nspin

        N_persite=count(h_io%onsite/=0.0d0)
        do i=1,N_neighbours
            N_persite=N_persite+count(h_io%hopping(:,:,i)/=0.0d0)*indexNN(i)
        enddo
        nnz=N_persite*Ncell

        !set the hamiltonian in coordinate matrix storage
        allocate(rowind(nnz),colind(nnz),source=0)
        allocate(val(nnz),source=cmplx(0.0d0,0.0d0,8))

        !onsite terms
        ii=0
        do i_orb=1,h_io%nb_orbitals
            ind_orb=(i_orb-1)*h_io%nb_spin
            do i_s=1,h_io%nb_spin
                if(h_io%onsite(i_s,i_orb)==0.0d0) cycle
                H_entry=cmplx(h_io%onsite(i_s,i_orb),0.0d0,8)
                do i_cell=1,Ncell
                    ii=ii+1
                    val(ii)=H_entry
                    rowind(ii)=i_s+ind_orb+(i_cell-1)*dim_mode
                    colind(ii)=rowind(ii)
                enddo
            enddo
        enddo

        !super ugly with tables, translate i_x,i_y,i_z to one index?
        !hopping terms
        do i_neigh=1,N_neighbours
            offset=sum(indexNN(1:i_neigh-1))
            do i_orb=1,h_io%nb_orbitals
                ind_orb=(i_orb-1)*h_io%nb_spin
                do i_s=1,h_io%nb_spin
                    if(h_io%hopping(i_s,i_orb,i_neigh)==0.0d0) cycle
                    H_entry=cmplx(h_io%hopping(i_s,i_orb,i_neigh),0.0d0,8)
                    do i_cell=1,Ncell
                        ilat_1=lat%index_1_3(i_cell)
                        do i_vois=1,indexNN(i_neigh)
                            if(tableNN(5,i_vois,ilat_1(1),ilat_1(2),ilat_1(3),1)/=1) cycle
                            ilat_2=tableNN(1:3,i_vois+offset,ilat_1(1),ilat_1(2),ilat_1(3),1)
                            i_other=lat%index_m_1(ilat_2)
                            ii=ii+1
                            val(ii)=H_entry
                            rowind(ii)=i_s+ind_orb+(i_cell-1)*dim_mode
                            colind(ii)=i_s+ind_orb+(i_other-1)*dim_mode
                        enddo
                    enddo
                enddo
            enddo
        enddo
        cols=dimH;rows=dimH
        stat = mkl_sparse_z_create_coo ( H_coo , SPARSE_INDEX_BASE_ONE , rows , cols , nnz , rowind , colind , val)
        if(stat /= 0) ERROR STOP "error creating H_coo in sparse H_ee"
        stat = MKL_SPARSE_CONVERT_CSR(H_coo,SPARSE_OPERATION_NON_TRANSPOSE,H_csr)
        if(stat /= 0) ERROR STOP "error setting H_ee sparse to CSR"
        stat=MKL_SPARSE_DESTROY(H_coo)
    end subroutine

    subroutine set_Hr_Jsd(h_par,mode_mag,Jsd,H_csr)
        !use m_eigen_interface, only : eigen_set_H_e_jsd
        !adds the Jsd coupling to a local real-space Hamiltonian Hr
        type(parameters_TB_Hsolve),intent(in)     ::  h_par
        type(SPARSE_MATRIX_T),intent(out)         ::  H_csr
        real(8),intent(in)           ::  Jsd(:)
        real(8),intent(in)                       ::  mode_mag(:,:)

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
        if(size(mode_mag,1)>3) ERROR STOP "will not work with nmag>1"

        ii=1
        do i=1,h_par%ncell
            do j=1,h_par%norb
                matind=[2*j-1,2*j]+(i-1)*h_par%norb*h_par%nspin
                jsdmat(1,1)=Jsd(j)*cmplx( mode_mag(3,i), 0.0d0        ,8)
                jsdmat(2,1)=Jsd(j)*cmplx( mode_mag(1,i), mode_mag(2,i),8)
                jsdmat(1,2)=Jsd(j)*cmplx( mode_mag(1,i),-mode_mag(2,i),8)
                jsdmat(2,2)=Jsd(j)*cmplx(-mode_mag(3,i), 0.0d0        ,8)
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
        if(stat /= 0) ERROR STOP "error creating H_coo for sparse H_Jsd"
        stat = MKL_SPARSE_CONVERT_CSR(H_coo,SPARSE_OPERATION_NON_TRANSPOSE,H_csr)
        if(stat /= 0) ERROR STOP "error setting H_Jsd sparse to CSR"
        stat=MKL_SPARSE_DESTROY(H_coo)
    end subroutine 
end module
#endif
