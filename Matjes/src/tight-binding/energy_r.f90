module m_energy_r
    use m_basic_types, only : vec_point
    use m_fermi,only: fermi_distrib
    use m_tb_types
    use m_energy_solve_dense
    use m_energy_set_real_sparse, only: set_Hr_sparse_nc
    use m_energy_set_real_sparse_sc, only: set_Hr_sparse_sc
    use m_energy_set_real_sc, only: set_Hr_dense_sc
    use m_energy_set_real, only: set_Hr_dense_nc
#ifdef CPP_MKL
    use MKL_SPBLAS
    use m_energy_solve_sparse
#endif

    implicit none

    private
    public :: get_eigenval_r,get_eigenvec_r,calc_occupation,set_Hr,set_Hr_dense
    !large electronic Hamiltonian
    complex(8),allocatable  ::  Hr(:,:)
#ifdef CPP_MKL
    type(SPARSE_MATRIX_T)   ::  Hr_sparse
    public :: set_hr_sparse
#endif

    contains

    subroutine calc_occupation(eigvec,eigval,E_f,kt,occ)
        real(8),intent(in)      ::   eigval(:),E_f,kt
        complex(8),intent(in)   ::   eigvec(:,:)
        real(8),intent(out)     ::   occ(size(eigval))
        integer                 ::  i

        occ=0.0d0
        do i=1,size(eigval)
            occ=occ+fermi_distrib(E_f,eigval(i),kt)*real(conjg(eigvec(:,i))*eigvec(:,i),kind=8)
        enddo
    end subroutine

    subroutine set_Hr(h_par,mode_mag)
        type(parameters_TB_Hsolve),intent(in)    ::  h_par
        type(vec_point),intent(in)               ::  mode_mag(:)

        if(h_par%sparse)then
#ifdef CPP_MKL
            Call set_Hr_sparse(h_par,mode_mag,Hr_sparse)
#else
            STOP "requires CPP_MKL for sparse TB"
#endif
        else
            Call set_Hr_dense(h_par,mode_mag,Hr)
        endif
    end subroutine 


#ifdef CPP_MKL
    subroutine set_Hr_sparse(h_par,mode_mag,Hr_set)
        type(parameters_TB_Hsolve),intent(in)    ::  h_par
        type(SPARSE_MATRIX_T),intent(out)        ::  Hr_set
        type(vec_point),intent(in)               ::  mode_mag(:)

        if(h_par%nsc==2)then
            Call set_Hr_sparse_sc(h_par,mode_mag,Hr_set)
        else
            Call set_Hr_sparse_nc(h_par,mode_mag,Hr_set)
        endif
    end subroutine 
#endif

    subroutine set_Hr_dense(h_par,mode_mag,Hr_set)
        type(parameters_TB_Hsolve),intent(in)    ::  h_par
        complex(8),allocatable,intent(inout)     ::  Hr_set(:,:)
        type(vec_point),intent(in)               ::  mode_mag(:)

        if(h_par%nsc==2)then
            Call set_Hr_dense_sc(h_par,mode_mag,Hr_set)
        else
            Call set_Hr_dense_nc(h_par,mode_mag,Hr_set)
        endif
    end subroutine 

    subroutine get_Hr(dimH,Hr_out)
        !not sure if needed anymore
        integer,intent(in)          ::  dimH
        complex(8),intent(out)      ::  Hr_out(dimH,dimH)
        
        if(.not.allocated(Hr)) STOP "Hr is not allocated but get_Hr is called"
        if(dimH/=size(Hr,1)) STOP "dimensions of Hr seems to be wrong for getting Hr"
        Hr_out=Hr
    end subroutine 

    subroutine get_eigenval_r(h_par,eigval,mode_mag)
        type(parameters_TB_Hsolve),intent(in)     ::  h_par
        real(8),intent(out),allocatable           ::  eigval(:)
        type(vec_point),intent(in)                ::  mode_mag(:)
    
        if(h_par%sparse)then
#ifdef CPP_MKL
            Call set_Hr_sparse(h_par,mode_mag,Hr_sparse)
            Call HR_eigval_sparse_feast(h_par,Hr_sparse,eigval)
#else
            STOP 'Cannot use spase get_eigenvalue without CPP_MKL'
#endif
        else
            Call set_Hr_dense(h_par,mode_mag,Hr)
            if(h_par%i_diag==1)then
                Call Hr_eigval_feast(h_par,Hr,eigval)
            elseif(h_par%i_diag==2)then
                Call Hr_eigval_zheev(h_par,Hr,eigval)
            elseif(h_par%i_diag==3)then
                Call Hr_eigval_zheevd(h_par,Hr,eigval)
            elseif(h_par%i_diag==4)then
                Call Hr_eigval_zheevr(h_par,Hr,eigval)
            else
                write(*,*) "trying to use h_par%i_diag=",h_par%i_diag
                STOP "h_par%diag choice not implemented for get_eigenval_r"
            endif
        endif
    end subroutine


    subroutine get_eigenvec_r(h_par,eigval,eigvec,mode_mag)
        type(parameters_TB_Hsolve),intent(in)     :: h_par
        real(8),intent(out),allocatable           :: eigval(:)
        complex(8),intent(out),allocatable        :: eigvec(:,:)
        type(vec_point),intent(in)                :: mode_mag(:)

        if(h_par%sparse)then
#ifdef CPP_MKL
            Call set_Hr_sparse(h_par,mode_mag,Hr_sparse)
            Call HR_eigvec_sparse_feast(h_par,Hr_sparse,eigvec,eigval)
#else
            STOP 'Cannot use spase get_eigenvalue without CPP_MKL'
#endif
        else
            Call set_Hr_dense(h_par,mode_mag,Hr)
            if(h_par%i_diag==1)then
                Call Hr_eigvec_feast(h_par,Hr,eigvec,eigval)
            elseif(h_par%i_diag==2)then
                Call Hr_eigvec_zheev(h_par,Hr,eigvec,eigval)
            elseif(h_par%i_diag==3)then
                Call Hr_eigvec_zheevd(h_par,Hr,eigvec,eigval)
            elseif(h_par%i_diag==4)then
                Call Hr_eigvec_zheevr(h_par,Hr,eigvec,eigval)
            else
                write(*,*) "trying to use h_par%i_diag=",h_par%i_diag
                STOP "h_par%diag choice not implemented for get_eigenval_r"
            endif
        endif

    end subroutine


end module m_energy_r
