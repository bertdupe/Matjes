module m_energy_r
    use m_basic_types, only : vec_point
    use m_fermi,only: fermi_distrib
    use m_tb_types
    use m_energy_solve_dense
!    use m_energy_set_real_sparse, only: set_Hr,set_Jsd,Hr_eigval,Hr_eigvec
    use m_energy_set_real_sc, only: set_Hr_dense_sc
    use m_energy_set_real, only: set_Hr_dense_nc

    implicit none

    !large electronic Hamiltonian
    complex(8),allocatable  ::  Hr(:,:)

    private
    public :: get_eigenval_r,get_eigenvec_r,calc_occupation,set_Hr
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

    subroutine set_Hr(h_par,mode_mag,Hr_set)
        type(parameters_TB_Hsolve),intent(in)     ::  h_par
        complex(8),allocatable,intent(inout)     ::  Hr_set(:,:)
        type(vec_point),intent(in)               ::  mode_mag(:)

        if(h_par%sparse)then
            STOP 'not implemented'
        else
            if(h_par%nsc==2)then
                Call set_Hr_dense_sc(h_par,mode_mag,Hr_set)
            else
                Call set_Hr_dense_nc(h_par,mode_mag,Hr_set)
            endif
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
        real(8),intent(out)         :: eigval(h_par%dimH)
        type(vec_point),intent(in)  :: mode_mag(:)

        Call set_Hr(h_par,mode_mag,Hr)
        !if(.false.)then
        !    Call Hr_eigval_feast(h_par%dimH,Hr,eigval)
        if(h_par%i_diag==2)then
            Call Hr_eigval_zheev(h_par%dimH,Hr,eigval)
        elseif(h_par%i_diag==3)then
            Call Hr_eigval_zheevd(h_par%dimH,Hr,eigval)
        else
            write(*,*) "trying to use h_par%i_diag=",h_par%i_diag
            STOP "h_par%diag choice not implemented for get_eigenval_r"
        endif
    end subroutine


    subroutine get_eigenvec_r(h_par,eigval,eigvec,mode_mag)
        type(parameters_TB_Hsolve),intent(in)     ::  h_par
        real(8),intent(out)         :: eigval(h_par%dimH)
        complex(8),intent(out)      :: eigvec(h_par%dimH,h_par%dimH)
        type(vec_point),intent(in)  :: mode_mag(:)

        Call set_Hr(h_par,mode_mag,Hr)
        if(h_par%i_diag==1)then
            Call Hr_eigvec_feast(h_par%dimH,Hr,eigvec,eigval)
        elseif(h_par%i_diag==2)then
            Call Hr_eigvec_zheev(h_par%dimH,Hr,eigvec,eigval)
        elseif(h_par%i_diag==3)then
            Call Hr_eigvec_zheevd(h_par%dimH,Hr,eigvec,eigval)
        else
            write(*,*) "trying to use h_par%i_diag=",h_par%i_diag
            STOP "h_par%diag choice not implemented for get_eigenval_r"
        endif

    end subroutine


end module m_energy_r
