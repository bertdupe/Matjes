module m_energy_r
    use m_basic_types, only : vec_point
    use m_fermi,only: fermi_distrib
    use m_rw_TB, only : TB_params
#ifdef CPP_SPARSE_HE__
    use m_energy_set_real_sparse, only: set_Hr,set_Jsd,Hr_eigval,Hr_eigvec
#elif CPP_SC
    use m_dos_sc, only : write_dos_sc 
    use m_energy_set_real_sc, only: set_Hr,Hr_eigval,Hr_eigvec
#else
    use m_energy_set_real, only: set_Hr,Hr_eigval,Hr_eigvec
#endif

    implicit none

    private
    public :: get_eigenval_r,get_occupation
    contains

    subroutine get_occupation(dimH,tb_ext,occ,mode_mag,E_f)
        integer,intent(in)          :: dimH
        integer,intent(in)          :: tb_ext(2)
        real(8),intent(out)         :: occ(dimH)
        real(8),intent(in)          :: E_f
        type(vec_point),intent(in)  :: mode_mag(dimH)

        complex(8)                  :: eigvec(dimH,dimH)
        real(8)                     :: eigval(dimH)
        
        Call set_Hr(dimH,tb_ext,mode_mag)
        Call Hr_eigvec(dimH,eigvec,eigval)
#if CPP_SC
        Call write_dos_sc(eigval,eigvec,'dos_r_sc.dat')
#endif
        Call calc_occupation(dimH,eigvec,eigval,E_f,TB_params%kt,occ)
    end subroutine 

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




end module m_energy_r
