module m_energy_r
    use m_basic_types, only : vec_point
    use m_fermi,only: fermi_distrib
    use m_rw_TB, only : TB_params
    use m_energy_set_real, only: set_Hr,H_add_Jsd

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

        real(8)                     :: eigval(dimH)
        complex(8)                  :: Hr(dimH,dimH)
        
        Call set_Hr(dimH,Hr,tb_ext)
        if(any(TB_params%Jsd /= 0.0d0))then
            Call H_add_Jsd(dimH,Hr,tb_ext,mode_mag,TB_params%Jsd)
        endif

        Call get_eigvec(dimH,Hr,eigval)
        Call calc_occupation(dimH,Hr,eigval,E_f,TB_params%kt,occ)
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

        complex(8)          :: Hr(dimH,dimH)
        
        Call set_Hr(dimH,Hr,tb_ext)
        if(any(TB_params%Jsd /= 0.0d0))then
            Call H_add_Jsd(dimH,Hr,tb_ext,mode_mag,TB_params%Jsd)
        endif
        Call get_eigval(dimH,Hr,eigval)
    end subroutine


    subroutine get_eigval(dimH,Hr,eigval)
        integer,intent(in)          ::  dimH
        complex(8),intent(inout)    ::  Hr(dimH,dimH)
        real(8),intent(out)         ::  eigval(dimH)

        complex(kind=8)             :: WORK(2*dimH)
        real(kind=8)                :: RWORK(3*dimH-2)
        integer                     :: info

        call ZHEEV( 'N', 'U', dimH, Hr, dimH, eigval, WORK, size(Work), RWORK, INFO )

    end subroutine

    subroutine get_eigvec(dimH,Hr,eigval)
        integer,intent(in)          ::  dimH
        complex(8),intent(inout)    ::  Hr(dimH,dimH)
        real(8),intent(out)         ::  eigval(dimH)

        complex(kind=8)             :: WORK(2*dimH)
        real(kind=8)                :: RWORK(3*dimH-2)
        integer                     :: info

        call ZHEEV( 'V', 'U', dimH, Hr, dimH, eigval, WORK, size(Work), RWORK, INFO )

    end subroutine 


end module m_energy_r
