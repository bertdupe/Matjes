module m_energy_r
    use m_basic_types, only : vec_point
    use m_fermi,only: fermi_distrib
    use m_tb_types
#ifdef CPP_SPARSE_HE__
    use m_energy_set_real_sparse, only: set_Hr,set_Jsd,Hr_eigval,Hr_eigvec
#elif CPP_SC
    use m_energy_set_real_sc, only: set_Hr,Hr_eigval,Hr_eigvec
#else
    use m_energy_set_real, only: set_Hr,Hr_eigval,Hr_eigvec
#endif

    implicit none

    private
    public :: get_eigenval_r,get_eigenvec_r,calc_occupation
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


    subroutine get_eigenval_r(Hsize,eigval,mode_mag)
        type(parameters_TB_Hsize),intent(in)     ::  Hsize
        real(8),intent(out)         :: eigval(Hsize%dimH)
        type(vec_point),intent(in)  :: mode_mag(:)

        Call set_Hr(Hsize,mode_mag)
        Call Hr_eigval(Hsize%dimH,eigval)
    end subroutine


    subroutine get_eigenvec_r(Hsize,eigval,eigvec,mode_mag)
        type(parameters_TB_Hsize),intent(in)     ::  Hsize
        real(8),intent(out)         :: eigval(Hsize%dimH)
        complex(8),intent(out)      :: eigvec(Hsize%dimH,Hsize%dimH)
        type(vec_point),intent(in)  :: mode_mag(:)

        Call set_Hr(Hsize,mode_mag)
        Call Hr_eigvec(Hsize%dimH,eigvec,eigval)
    end subroutine


end module m_energy_r
