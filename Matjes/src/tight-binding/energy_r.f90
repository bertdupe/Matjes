module m_energy_r
    use m_basic_types, only : vec_point
    use m_fermi,only: fermi_distrib
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

    subroutine calc_occupation(dimH,eigvec,eigval,E_f,kt,occ)
        integer,intent(in)      ::   dimH
        real(8),intent(in)      ::   eigval(dimH),E_f,kt
        complex(8),intent(in)   ::   eigvec(dimH,dimH)
        real(8),intent(out)     ::   occ(dimH)

        integer                 ::  i

        occ=0.0d0
        do i=1,dimH
            occ=occ+fermi_distrib(E_f,eigval(i),kt)*real(conjg(eigvec(:,i))*eigvec(:,i),kind=8)
        enddo
    end subroutine


    subroutine get_eigenval_r(dimH,tb_ext,eigval,mode_mag)
        integer,intent(in)          :: dimH
        integer,intent(in)          :: tb_ext(2)
        real(8),intent(out)         :: eigval(dimH)
        type(vec_point),intent(in)  :: mode_mag(:)

        Call set_Hr(dimH,tb_ext,mode_mag)
        Call Hr_eigval(dimH,eigval)
    end subroutine


    subroutine get_eigenvec_r(dimH,tb_ext,eigval,eigvec,mode_mag)
        integer,intent(in)          :: dimH
        integer,intent(in)          :: tb_ext(2)
        real(8),intent(out)         :: eigval(dimH)
        complex(8),intent(out)      :: eigvec(dimH,dimH)
        type(vec_point),intent(in)  :: mode_mag(:)

        Call set_Hr(dimH,tb_ext,mode_mag)
        Call Hr_eigvec(dimH,eigvec,eigval)
    end subroutine


end module m_energy_r
