module m_H_tb_public
use m_H_tb_base
use m_H_tb_coo
use m_H_tb_dense
use m_H_tb_csr
use m_TB_types, only: parameters_TB_IO_H
private set_H_single, set_H_multiple
interface set_H
    module procedure set_H_single
    module procedure set_H_multiple
end interface

contains

subroutine restict_solution_positive(eval,evec,emax,vanish)
    !restrict the eigenvalues and eigenvectors only to positive energies, up to the optional maximal energy emax
    complex(8),allocatable,intent(inout)    :: evec(:,:) !eigenvectors
    real(8),allocatable,intent(inout)       :: eval(:)   !eigenvalues
    real(8),optional,intent(in)             :: emax
    logical,intent(out),optional            :: vanish 
    complex(8),allocatable  :: tmp_c(:,:)
    real(8),allocatable     :: tmp_r(:)

    integer     :: ie_bnd(2)
    integer     :: i

    if(present(vanish)) vanish=.false.

    !only sum over positive energies
    ie_bnd=size(eval)+1
    do i=1,size(eval)
        if(eval(i)>0.0d0)then
            ie_bnd(1)=i
            exit
        endif
    enddo
    if(ie_bnd(1)>size(eval)) STOP "need positive eigenvalues for selfconsistent delta"
    ie_bnd(2)=size(eval)
    if(present(emax).and.emax>0.0d0)then
        do i=ie_bnd(1),size(eval)
            if(eval(i)>emax)then
                ie_bnd(2)=i-1
                exit
            endif
        enddo
    endif
    if(ie_bnd(2)<ie_bnd(1))then
        if(present(vanish))then
            vanish=.true.
            deallocate(evec,eval)
            return
        else
            STOP "there needs to be a non-vanishing set of positive eigenvalues up to emax" !reference input parameter setting emax once defined
        endif
    endif
    Call move_alloc(evec,tmp_c)
    Call move_alloc(eval,tmp_r)
    allocate(eval,source=tmp_r(  ie_bnd(1):ie_bnd(2)))
    deallocate(tmp_r)
    allocate(evec,source=tmp_c(:,ie_bnd(1):ie_bnd(2)))
    deallocate(tmp_c)
end subroutine


subroutine H_append(H,Hadd)
    !adds the entries from Hadd to H, and destroys the Hadd array
    type(H_tb_coo),allocatable,intent(inout)   :: H(:)
    type(H_tb_coo),allocatable,intent(inout)   :: Hadd(:)

    type(H_tb_coo),allocatable                 :: Htmp(:)
    !unfortunately does not work with classes...
    !class(H_tb),allocatable,intent(inout)   :: H(:)
    !class(H_tb),allocatable,intent(inout)   :: Hadd(:)

    !class(H_tb),allocatable                 :: Htmp(:)

    integer ::  i

    if(.not.allocated(Hadd)) ERROR STOP "Cannot append H as the added H is not allocated"
    if(allocated(H))then
        allocate(Htmp(size(H)+size(Hadd)),mold=H)
        do i=1,size(H)
            Call H(i)%mv(Htmp(i))
        enddo
        do i=1,size(Hadd)
            Call Hadd(i)%mv(Htmp(i+size(H)))
        enddo
        deallocate(H,Hadd)
        Call move_alloc(Htmp,H)
    else
        Call move_alloc(Hadd,H)
    endif
end subroutine

subroutine set_H_single(H,io)
    use, intrinsic :: iso_fortran_env, only : output_unit, error_unit
    class(H_tb),allocatable,intent(inout)   :: H
    type(parameters_TB_IO_H),intent(in)     :: io
    logical,save    :: nsaid=.true.

    if(allocated(H)) STOP "CANNOT set H which is already set"
    select case(io%i_diag)
    case(1)
#ifdef CPP_LAPACK
        if(nsaid) write(output_unit,'(2/A/)') "Chose lapack zheevd algoritm for tight-binding Hamiltonian"
        allocate(H_zheevd::H)
#else
        write(error_unit,'(//A)') "CANNOT use lapack diagonalization algorithm without CPP_LAPACK"
        ERROR STOP
#endif
    case(2)
#ifdef CPP_LAPACK
            if(nsaid) write(output_unit,'(2/A/)') "Chose lapack zheev algoritm for tight-binding Hamiltonian"
            allocate(H_zheev::H)
#else
            write(error_unit,'(//A)') "CANNOT use lapack diagonalization algorithm without CPP_LAPACK"
            ERROR STOP
#endif
    case(3)
#ifdef CPP_MKL
        if(nsaid) write(output_unit,'(2/A/)') "Chose dense feast algoritm for tight-binding Hamiltonian"
        allocate(H_feast_den::H)
#else
        write(error_unit,'(//A)') "CANNOT use dense feast algorithm without CPP_MKL"
        write(error_unit,'(A)')   "Choose different TB_diag  or compile with USE_MKL"
        ERROR STOP 
#endif
    case(4)
#ifdef CPP_LAPACK
        if(nsaid) write(output_unit,'(2/A/)') "Chose lapack zheevr algoritm for tight-binding Hamiltonian"
        allocate(H_zheevr::H)
#else
        write(error_unit,'(//A)') "CANNOT use lapack diagonalization algorithm without CPP_LAPACK"
        ERROR STOP
#endif
    case(5)
#ifdef CPP_LAPACK
        if(nsaid) write(output_unit,'(2/A/)') "Chose lapack zheevx algoritm for tight-binding Hamiltonian"
        allocate(H_zheevx::H)
#else
        write(error_unit,'(//A)') "CANNOT use lapack diagonalization algorithm without CPP_LAPACK"
        ERROR STOP
#endif
    case(6)
#ifdef CPP_MKL
        write(output_unit,'(2/A/)') "Chose sparse feast algoritm for tight-binding Hamiltonian"
        allocate(H_feast_csr::H)
#else
        write(error_unit,'(//A)') "CANNOT use sparse feast algorithm without CPP_MKL"
        write(error_unit,'(A)')   "Disable sparse TB desription or compile with USE_MKL"
        ERROR STOP 
#endif
    case default
        write(error_unit,'(2/A,I6,A)') "Unable to choose dense tight-binding Hamiltonian as TB_diag=",io%i_diag," is not implemented"
        STOP "CHECK INPUT"
    end select
    nsaid=.false.
end subroutine

subroutine set_H_multiple(H,N,io)
    use, intrinsic :: iso_fortran_env, only : output_unit, error_unit
    class(H_tb),allocatable,intent(inout)   :: H(:)
    integer,intent(in)                      :: N
    type(parameters_TB_IO_H),intent(in)     :: io

    if(allocated(H)) STOP "CANNOT set H which is already set"
    select case(io%i_diag)
    case(1)
#ifdef CPP_LAPACK
        write(output_unit,'(2/A/)') "Chose lapack zheevd algoritm for tight-binding Hamiltonian"
        allocate(H_zheevd::H(N))
#else
        write(error_unit,'(//A)') "CANNOT use lapack diagonalization algorithm without CPP_LAPACK"
        ERROR STOP
#endif
    case(2)
#ifdef CPP_LAPACK
        write(output_unit,'(2/A/)') "Chose lapack zheev algoritm for tight-binding Hamiltonian"
        allocate(H_zheev::H(N))
#else
        write(error_unit,'(//A)') "CANNOT use lapack diagonalization algorithm without CPP_LAPACK"
        ERROR STOP
#endif
    case(3)
#ifdef CPP_MKL
        write(output_unit,'(2/A/)') "Chose dense feast algoritm for tight-binding Hamiltonian"
        allocate(H_feast_den::H(N))
#else
        write(error_unit,'(//A)') "CANNOT use dense feast algorithm without CPP_MKL"
        write(error_unit,'(A)')   "Choose different TB_diag  or compile with USE_MKL"
        ERROR STOP 
#endif
    case(4)
#ifdef CPP_LAPACK
        write(output_unit,'(2/A/)') "Chose lapack zheevr algoritm for tight-binding Hamiltonian"
        allocate(H_zheevr::H(N))
#else
        write(error_unit,'(//A)') "CANNOT use lapack diagonalization algorithm without CPP_LAPACK"
        ERROR STOP
#endif
    case(5)
#ifdef CPP_LAPACK
        write(output_unit,'(2/A/)') "Chose lapack zheevx algoritm for tight-binding Hamiltonian"
        allocate(H_zheevx::H(N))
#else
        write(error_unit,'(//A)') "CANNOT use lapack diagonalization algorithm without CPP_LAPACK"
        ERROR STOP
#endif
    case(6)
#ifdef CPP_MKL
        write(output_unit,'(2/A/)') "Chose sparse feast algoritm for tight-binding Hamiltonian"
        allocate(H_feast_csr::H(N))
#else
        write(error_unit,'(//A)') "CANNOT use sparse feast algorithm without CPP_MKL"
        write(error_unit,'(A)')   "Disable sparse TB desription or compile with USE_MKL"
        ERROR STOP 
#endif
    case default
        write(error_unit,'(2/A,I6,A)') "Unable to choose dense tight-binding Hamiltonian as TB_diag=",io%i_diag," is not implemented"
        STOP "CHECK INPUT"
    end select
end subroutine


end module
