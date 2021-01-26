module m_init_Hk
use m_derived_types, only: lattice
use m_H_tb_public
use m_tb_types ,only: parameters_TB_Hsolve,parameters_TB_IO_H, parameters_ham_init 
use m_types_tb_h_inp 
use m_neighbor_type, only: neighbors
use m_init_H
implicit none

private
public get_Hk_inp, Hk_inp_t, Hk_evec, HK_eval

type Hk_inp_t
    type(H_tb_coo),allocatable  :: H(:)
    real(8),allocatable         :: diffR(:,:)
end type

contains
subroutine get_Hk_inp(lat,h_io,Hk_inp)
    use m_init_HR, only: get_delta
    !get Hk_inp, which contains the Hamiltonian data to construct the Hamiltonian at arbitrary k
    type(lattice),intent(in)                :: lat
    type(parameters_TB_IO_H),intent(in)     :: h_io
    type(Hk_inp_t),intent(out)              :: Hk_inp

    type(parameters_TB_IO_H)                :: h_io_scf

    if(h_io%use_scf)then    
        if(.true.)then 
            !use reciprocal-space self-consistent delta

        else
            !use real-space self-consistent delta
            H_io_scf=h_io
            Call get_delta(lat,h_io,h_io_scf%del)
            Call get_H(lat,h_io_scf,Hk_inp%H,diffR=Hk_inp%diffR)
        endif
    else
        Call get_H(lat,h_io,Hk_inp%H,diffR=Hk_inp%diffR)
    endif
end subroutine

subroutine Hk_eval(Hk_inp,k,h_io,eval)
    type(Hk_inp_t),intent(in)               :: Hk_inp
    real(8),intent(in)                      :: k(3)
    type(parameters_TB_IO_H),intent(in)     :: h_io
    real(8),intent(out),allocatable         :: eval(:)
    class(H_tb),allocatable                 :: H

    Call get_Hk(Hk_inp,k,h_io,H)
    Call H%get_eval(eval)
    Call H%destroy
end subroutine

subroutine Hk_evec(Hk_inp,k,h_io,eval,evec)
    type(Hk_inp_t),intent(in)               :: Hk_inp
    real(8),intent(in)                      :: k(3)
    type(parameters_TB_IO_H),intent(in)     :: h_io
    real(8),allocatable,intent(out)         :: eval(:)
    complex(8),allocatable,intent(out)      :: evec(:,:)
    class(H_tb),allocatable                 :: H

    Call get_Hk(Hk_inp,k,h_io,H)
    Call H%get_evec(eval,evec)
    Call H%destroy
end subroutine

subroutine get_Hk(Hk_inp,k,h_io,H)
    !unfolds the Hk_inp data to and Hamiltonian (H) at a k-point (k)
    type(Hk_inp_t),intent(in)   :: Hk_inp
    real(8),intent(in)          :: k(3)
    type(parameters_TB_IO_H),intent(in)     :: h_io
    class(H_tb),allocatable,intent(out)     :: H
    class(H_tb),allocatable                 :: Htmp

    type(parameters_ham_init)   :: hinit
    integer     :: i_ham,N_ham
    complex(8),allocatable  ::  val(:)
    integer,allocatable     ::  row(:),col(:)

    real(8)  :: phase

    Call set_H(H,h_io)
    allocate(Htmp,mold=H)
    N_ham=size(HK_inp%H)
    do i_ham=1,N_ham
        phase=dot_product(HK_inp%diffR(:,i_ham),k)
        Call Hk_inp%H(i_ham)%get_hinit(hinit)
        Call Hk_inp%H(i_ham)%get_par(val,row,col)
        val=val*cmplx(cos(phase),sin(phase),8)
        Call Htmp%init_coo(val,row,col,hinit)
        Call H%add(Htmp)
        Call Htmp%destroy()
        deallocate(val,row,col)
    enddo
end subroutine
end module


